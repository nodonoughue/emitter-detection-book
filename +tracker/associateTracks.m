function [assigned_track_idx, assigned_msmt_idx, unassigned_msmt_idx, pda_states] = ...
    associateTracks(tracks, measurements, curr_time, motion_model, msmt_model, ...
                    gate_probability, assoc_type, detection_probability)
% associateTracks  Associate measurements with tracks (NN, GNN, or PDA).
%
% [assigned_track_idx, assigned_msmt_idx, unassigned_msmt_idx] = ...
%     associateTracks(tracks, measurements, curr_time, motion_model, msmt_model, ...
%                     gate_probability, assoc_type)
%
% [assigned_track_idx, assigned_msmt_idx, unassigned_msmt_idx, pda_states] = ...
%     associateTracks(..., detection_probability)
%
% INPUTS
%   tracks               Cell array of Track structs
%   measurements         Cell array of Measurement structs.  Each measurement's
%                        .msmt_model field is used; if empty, the explicit
%                        msmt_model argument is used as a fallback.
%   curr_time            Current scan timestamp [s] (used when measurements is empty)
%   motion_model         Motion model struct from makeMotionModel
%   msmt_model           Measurement model struct (fallback when msmt.msmt_model is empty)
%   gate_probability     Chi-square acceptance gate probability (e.g., 0.99)
%   assoc_type           'nn'  – Nearest Neighbour
%                        'gnn' – Global Nearest Neighbour (Munkres/Hungarian)
%                        'pda' – Probabilistic Data Association
%                        Default: 'gnn'
%   detection_probability  Pd, probability the target generates a measurement.
%                          Used only by PDA. Default: 1.0
%
% OUTPUTS
%   assigned_track_idx  (K x 1) index into tracks for each assignment
%   assigned_msmt_idx   (K x 1) index into measurements (0 = coast / PDA-handled)
%   unassigned_msmt_idx Indices of measurements consumed by no track (available
%                       to the initiator)
%   pda_states          (NN/GNN) empty cell {}.
%                       (PDA)    Cell array, one fused State struct per track.
%                                Callers should apply these directly rather than
%                                calling ekfUpdate.
%
% Nicholas O'Donoughue
% June 2025

if nargin < 7 || isempty(assoc_type)
    assoc_type = 'gnn';
end
if nargin < 8 || isempty(detection_probability)
    detection_probability = 1.0;
end

pda_states = {};   % empty for NN / GNN

num_tracks = numel(tracks);
num_msmts  = numel(measurements);

% If no measurements, all tracks coast
if num_msmts == 0
    assigned_track_idx  = (1:num_tracks)';
    assigned_msmt_idx   = zeros(num_tracks, 1);
    unassigned_msmt_idx = zeros(0, 1);
    return;
end

% Attach msmt_model fallback to any measurement that lacks one
for jj = 1:num_msmts
    if isempty(measurements{jj}.msmt_model)
        measurements{jj}.msmt_model = msmt_model;
    end
end

% Chi-square gate threshold (degrees of freedom = measurement dimension)
num_msmt_dims = numel(measurements{1}.zeta);
gate_size     = chi2inv(gate_probability, num_msmt_dims);

% Solve the assignment problem
switch lower(assoc_type)
    case 'nn'
        dist_mat = build_dist_mat(tracks, measurements, motion_model, gate_size, ...
                                  num_tracks, num_msmts);
        [assigned_track_idx, assigned_msmt_idx] = nn_assign(dist_mat);
        unassigned_msmt_idx = setdiff(1:num_msmts, assigned_msmt_idx(assigned_msmt_idx > 0))';

    case 'gnn'
        dist_mat = build_dist_mat(tracks, measurements, motion_model, gate_size, ...
                                  num_tracks, num_msmts);
        [assigned_track_idx, assigned_msmt_idx] = gnn_assign(dist_mat, gate_size);
        unassigned_msmt_idx = setdiff(1:num_msmts, assigned_msmt_idx(assigned_msmt_idx > 0))';

    case 'pda'
        [pda_states, unassigned_msmt_idx] = pda_assign( ...
            tracks, measurements, curr_time, motion_model, ...
            gate_size, gate_probability, detection_probability);
        assigned_track_idx  = (1:num_tracks)';
        % msmt_idx = 0 for all tracks: the caller must use pda_states instead
        % of running a standard ekfUpdate.
        assigned_msmt_idx   = zeros(num_tracks, 1);

    otherwise
        error('tracker:associateTracks:unknownType', ...
              'Unknown association type "%s"; use ''nn'', ''gnn'', or ''pda''.', assoc_type);
end

end


%% ---- Shared distance-matrix builder (NN and GNN) -------------------------

function dist_mat = build_dist_mat(tracks, measurements, motion_model, gate_size, ...
                                   num_tracks, num_msmts)
% Build (num_tracks x num_msmts) distance matrix.
% Entries outside the chi-square gate are set to Inf.
dist_mat = inf(num_tracks, num_msmts);
for ii = 1:num_tracks
    trk = tracks{ii};
    if isempty(trk.motion_model)
        trk.motion_model = motion_model;
    end
    for jj = 1:num_msmts
        h = tracker.makeHypothesis(trk, measurements{jj});
        if h.distance <= gate_size
            dist_mat(ii, jj) = h.distance;
        end
    end
end
end


%% ---- NN assignment --------------------------------------------------------

function [trk_idx, msmt_idx] = nn_assign(dist_mat)
% Sequential nearest-neighbour assignment.
num_tracks = size(dist_mat, 1);
trk_idx  = (1:num_tracks)';
msmt_idx = zeros(num_tracks, 1);
used     = false(1, size(dist_mat, 2));

for ii = 1:num_tracks
    row = dist_mat(ii, :);
    row(used) = inf;
    [d_min, jj] = min(row);
    if isfinite(d_min)
        msmt_idx(ii) = jj;
        used(jj) = true;
    end
end
end


%% ---- GNN assignment -------------------------------------------------------

function [trk_idx, msmt_idx] = gnn_assign(dist_mat, null_cost)
% Global nearest-neighbour assignment via the Hungarian (Munkres) algorithm.
% Augment cost matrix with one null column per track at cost null_cost
% (= chi2inv gate threshold), so a track prefers a null column over any
% out-of-gate measurement (cost = big >> null_cost).
num_tracks = size(dist_mat, 1);
num_msmts  = size(dist_mat, 2);

% Replace inf (out-of-gate) with a cost much larger than null_cost.
big  = null_cost * 1e6;
cost = dist_mat;
cost(~isfinite(cost)) = big;

% Pad with null columns (one per track) at the gate threshold cost.
cost_aug = [cost, null_cost * ones(num_tracks, num_tracks)];

if exist('assignDetectionsToTracks', 'file') == 2
    % Automated Driving Toolbox: non-assignment cost threshold = null_cost.
    [assignments, ~, ~] = assignDetectionsToTracks(cost, null_cost);
    trk_idx  = (1:num_tracks)';
    msmt_idx = zeros(num_tracks, 1);
    for k = 1:size(assignments, 1)
        msmt_idx(assignments(k, 1)) = assignments(k, 2);
    end
else
    % Fallback: munkres with augmented null columns.
    [row_assign, ~] = munkres(cost_aug);
    trk_idx  = (1:num_tracks)';
    msmt_idx = zeros(num_tracks, 1);
    for ii = 1:num_tracks
        jj = row_assign(ii);
        if jj <= num_msmts && isfinite(dist_mat(ii, jj))
            msmt_idx(ii) = jj;   % real in-gate measurement
        end
        % jj > num_msmts, or originally Inf (out-of-gate) → coast
    end
end
end


%% ---- PDA assignment -------------------------------------------------------

function [pda_states, ungated_msmt_idx] = pda_assign( ...
        tracks, measurements, curr_time, motion_model, ...
        gate_size, gate_probability, detection_probability)
% Probabilistic Data Association update.
%
% For each track, all gated measurements contribute to the EKF update via
% normalised Gaussian-likelihood weights.  A null hypothesis (no detection)
% is always included with weight p_miss = 1 - Pd * Pg.  The returned state
% is the Gaussian mixture mean with mixture covariance (including the
% spread-of-means term).
%
% This is equivalent to Python's PDAAssociator.associate() followed by
% GMMHypothesis.update_track() on every track.
%
% INPUTS
%   tracks             Cell array of Track structs
%   measurements       Cell array of Measurement structs (all same scan time)
%   curr_time          Timestamp to predict to [s]
%   motion_model       Motion model struct
%   gate_size          Chi-square gate threshold (chi2inv(gate_probability, n_msmt_dim))
%   gate_probability   Pg — used to compute p_miss = 1 - Pd * Pg
%   detection_probability  Pd (default 1.0)
%
% OUTPUTS
%   pda_states        Cell array, one fused State struct per track
%   ungated_msmt_idx  Row vector of measurement indices not gated by any track

num_tracks    = numel(tracks);
num_msmts     = numel(measurements);
n_msmt_dim    = numel(measurements{1}.zeta);
p_miss        = 1 - detection_probability * gate_probability;

pda_states   = cell(num_tracks, 1);
gated_by_any = false(1, num_msmts);

for ii = 1:num_tracks
    % Measurement model (shared across all measurements in the scan)
    mm = measurements{1}.msmt_model;

    % Resolve motion model: prefer track-carried model, fall back to argument
    trk_mm = tracks{ii}.motion_model;
    if isempty(trk_mm)
        trk_mm = motion_model;
    end

    % Predict track to current scan time
    s      = tracker.currState(tracks{ii});
    s_pred = tracker.predictState(s, curr_time, trk_mm);

    % Linearize at the predicted state — shared for all measurements
    z_hat = mm.z_fun(s_pred);
    H     = mm.h_fun(s_pred);
    if ndims(H) == 3
        H = H(:, :, 1)';    % squeeze (n_st x n_m x 1) → (n_m x n_st)
    end
    P = s_pred.covar;
    R = mm.R;
    S = H * P * H' + R;
    K = P * H' / S;
    n_st = numel(s_pred.state);

    % Gate measurements and collect Gaussian likelihoods
    gated_idx  = [];
    gated_like = [];
    for jj = 1:num_msmts
        nu = measurements{jj}.zeta(:) - z_hat(:);
        d2 = nu' / S * nu;
        if d2 <= gate_size
            % N(nu ; 0, S) evaluated analytically (no Statistics Toolbox needed)
            L_j = exp(-0.5 * d2) / sqrt((2*pi)^n_msmt_dim * det(S));
            gated_idx(end+1)  = jj;   %#ok<AGROW>
            gated_like(end+1) = L_j;  %#ok<AGROW>
            gated_by_any(jj) = true;
        end
    end

    % Normalised weights: [L_1 ... L_M  p_miss] / total
    all_weights = [gated_like, p_miss];
    all_weights = all_weights / sum(all_weights);

    % Per-hypothesis state updates
    n_gated = numel(gated_idx);
    n_hyp   = n_gated + 1;
    all_states = zeros(n_st, n_hyp);
    all_covars = zeros(n_st, n_st, n_hyp);

    I_nst = eye(n_st);
    for kk = 1:n_gated
        nu_j       = measurements{gated_idx(kk)}.zeta(:) - z_hat(:);
        all_states(:, kk) = s_pred.state + K * nu_j;
        P_upd = (I_nst - K * H) * P;
        % Enforce symmetry and PSD (matches Python CovarianceMatrix.ensure_positive_definite)
        P_upd = (P_upd + P_upd') / 2;
        [V_u, D_u] = eig(P_upd);
        P_upd = V_u * diag(max(diag(D_u), 1e-10)) * V_u';
        all_covars(:, :, kk) = (P_upd + P_upd') / 2;
    end
    % Null hypothesis: coast (no measurement update)
    all_states(:, end)    = s_pred.state;
    all_covars(:, :, end) = P;

    % Fused mean
    x_fused = all_states * all_weights(:);

    % Fused covariance: weighted sum of (per-hypothesis covar + spread-of-means)
    P_fused = zeros(n_st);
    for kk = 1:n_hyp
        delta   = all_states(:, kk) - x_fused;
        P_fused = P_fused + all_weights(kk) * (all_covars(:, :, kk) + delta * delta');
    end

    pda_states{ii} = tracker.makeState(s_pred.state_space, s_pred.time, x_fused, P_fused);
end

% Measurements not gated by any track are available to the initiator
ungated_msmt_idx = find(~gated_by_any)';
end


%% ---- Munkres (Hungarian) algorithm ---------------------------------------

function [row_ind, col_ind] = munkres(costMatrix)
% Munkres (Hungarian) algorithm.
%
% Structured to match the standard four-step description:
%   Step 2 – cover starred columns; done if enough covered.
%   Step 3 – scan for uncovered unstarred zeros; prime each one found.
%             If no star in that row  → Step 5 (augmenting path).
%             If star in that row     → cover row, uncover star's column,
%                                       repeat Step 3.
%             If no uncovered zero    → Step 6 (cost adjust), then Step 3.
%   Step 5 – collect augmenting path and flip; reset covers → Step 2.
%   Step 6 – subtract delta from uncovered, add to doubly-covered → Step 3.
%
% The key correctness point: Step 6 stays in the inner loop so that
% covered_rows / covered_cols are preserved across cost adjustments.
% Only Step 5 breaks back to Step 2 and resets the covers.

[n, m] = size(costMatrix);
C = costMatrix;
C = bsxfun(@minus, C, min(C, [], 2));   % subtract row minima
C = bsxfun(@minus, C, min(C, [], 1));   % subtract col minima

star  = false(n, m);
prime = false(n, m);

% Initial greedy starring
for i = 1:n
    for j = 1:m
        if C(i,j) == 0 && ~any(star(i,:)) && ~any(star(:,j))
            star(i,j) = true;
        end
    end
end

while true
    %% Step 2: cover columns that contain a starred zero
    covered_cols = any(star, 1);            % reset every time we reach Step 2
    covered_rows = false(n, 1);
    if sum(covered_cols) >= min(n, m)
        break;                              % optimal assignment found
    end

    %% Steps 3 / 5 / 6 inner loop (covers persist across Step 6)
    while true
        %% Step 3: find an uncovered, unstarred zero (vectorised)
        zero_mask = (C == 0) & ~star & ~covered_rows & ~covered_cols;
        [z_row, z_col] = find(zero_mask, 1, 'first');
        found_zero = ~isempty(z_row);
        if ~found_zero, z_row = 0; z_col = 0; end

        if ~found_zero
            %% Step 6: cost adjustment — stay in inner loop
            uncov = C(~covered_rows, ~covered_cols);
            delta = min(uncov(:));
            C(~covered_rows, :) = C(~covered_rows, :) - delta;
            C(:, covered_cols)  = C(:, covered_cols)  + delta;
            continue;   % back to Step 3 with covers intact
        end

        prime(z_row, z_col) = true;
        sc = find(star(z_row, :), 1);

        if ~isempty(sc)
            %% Step 3 continued: star in this row — cover row, uncover col
            covered_rows(z_row) = true;
            covered_cols(sc)    = false;
            continue;   % back to Step 3
        end

        %% Step 5: augmenting path from (z_row, z_col)
        curr_col = z_col;
        path_r   = z_row;
        path_c   = z_col;
        while true
            si = find(star(:, curr_col), 1);
            if isempty(si), break; end
            pj2      = find(prime(si, :), 1);
            path_r(end+1) = si;   path_c(end+1) = curr_col;  %#ok<AGROW>
            path_r(end+1) = si;   path_c(end+1) = pj2;       %#ok<AGROW>
            curr_col = pj2;
        end
        for k = 1:2:numel(path_r)
            star(path_r(k), path_c(k)) = true;
        end
        for k = 2:2:numel(path_r)
            star(path_r(k), path_c(k)) = false;
        end
        prime(:) = false;
        break;  % back to Step 2 (outer loop resets covers)
    end
end

% Extract row → col assignment
row_ind = zeros(1, n);
for i = 1:n
    j = find(star(i,:), 1);
    if ~isempty(j)
        row_ind(i) = j;
    end
end
col_ind = row_ind;
end
