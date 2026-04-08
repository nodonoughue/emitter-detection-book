function [assigned_track_idx, assigned_msmt_idx, unassigned_msmt_idx] = ...
    associateTracks(tracks, measurements, curr_time, motion_model, msmt_model, gate_probability, assoc_type)
% associateTracks  Associate measurements with tracks (NN or GNN).
%
% [assigned_track_idx, assigned_msmt_idx, unassigned_msmt_idx] = ...
%     associateTracks(tracks, measurements, curr_time, motion_model, msmt_model, gate_probability, assoc_type)
%
% INPUTS
%   tracks          Cell array of Track structs
%   measurements    Cell array of Measurement structs
%   curr_time       Current scan timestamp [s] (used for coasting when measurements is empty)
%   motion_model    Motion model struct from makeMotionModel
%   msmt_model      Measurement model struct from makeMeasurementModel
%   gate_probability Chi-square acceptance gate probability (e.g., 0.99)
%   assoc_type      'nn'  – Nearest Neighbour (sequential, priority to first tracks)
%                   'gnn' – Global Nearest Neighbour (Munkres/Hungarian, globally optimal)
%                   Default: 'gnn'
%
% OUTPUTS
%   assigned_track_idx  (K x 1) integer index into tracks for each assignment
%   assigned_msmt_idx   (K x 1) integer index into measurements for each assignment
%                       0 means the track received a coast (missed detection)
%   unassigned_msmt_idx (M-K x 1) indices of measurements not assigned to any track
%
% Nicholas O'Donoughue
% June 2025

if nargin < 7 || isempty(assoc_type)
    assoc_type = 'gnn';
end

num_tracks = numel(tracks);
num_msmts  = numel(measurements);

% Pre-compute predicted states and gate sizes
num_msmt_dims = numel(measurements{1}.zeta);
gate_size = chi2inv(gate_probability, num_msmt_dims);

% Build distance matrix (num_tracks x num_msmts)
dist_mat = inf(num_tracks, num_msmts);
for ii = 1:num_tracks
    s = tracker.currState(tracks{ii});
    s_pred = tracker.predictState(s, curr_time, motion_model);
    for jj = 1:num_msmts
        d = tracker.computeDistance(s_pred, measurements{jj}.zeta, msmt_model);
        if d * num_msmt_dims <= gate_size  % un-normalise to compare with chi2 gate
            dist_mat(ii, jj) = d;
        end
    end
end

% Solve the assignment problem
switch lower(assoc_type)
    case 'nn'
        [assigned_track_idx, assigned_msmt_idx] = nn_assign(dist_mat);
    case 'gnn'
        [assigned_track_idx, assigned_msmt_idx] = gnn_assign(dist_mat);
    otherwise
        error('tracker:associateTracks:unknownType', ...
              'Unknown association type "%s"; use ''nn'' or ''gnn''.', assoc_type);
end

% Tracks with no finite assignment get coast (index 0)
all_assigned_msmt = assigned_msmt_idx(assigned_msmt_idx > 0);
unassigned_msmt_idx = setdiff(1:num_msmts, all_assigned_msmt)';

end


%% ---- Local helpers -------------------------------------------------------

function [trk_idx, msmt_idx] = nn_assign(dist_mat)
% Sequential nearest-neighbour assignment.
num_tracks = size(dist_mat, 1);
trk_idx  = (1:num_tracks)';
msmt_idx = zeros(num_tracks, 1);  % 0 = no assignment (coast)
used = false(1, size(dist_mat, 2));

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


function [trk_idx, msmt_idx] = gnn_assign(dist_mat)
% Global nearest-neighbour assignment via the Hungarian (Munkres) algorithm.
% Augment cost matrix so unassigned tracks do not block each other.
num_tracks = size(dist_mat, 1);
num_msmts  = size(dist_mat, 2);

% Replace inf with a large cost; add dummy columns for missed detections
big = 1e9;
cost = dist_mat;
cost(~isfinite(cost)) = big;

% Pad with dummy measurement columns so every track can be assigned
cost_aug = [cost, big * ones(num_tracks, num_tracks)];

% Run Hungarian algorithm (MATLAB built-in assignmentoptimal or munkres)
% Use MATLAB's built-in assignDetectionsToTracks if the Automated Driving
% Toolbox is available; otherwise fall back to a simple implementation.
if exist('assignDetectionsToTracks', 'file') == 2
    % Automated Driving Toolbox available
    [assignments, unassigned_trk, ~] = assignDetectionsToTracks(cost_aug, big - 1);
    trk_idx  = (1:num_tracks)';
    msmt_idx = zeros(num_tracks, 1);
    for k = 1:size(assignments, 1)
        t = assignments(k, 1);
        m = assignments(k, 2);
        if m <= num_msmts
            msmt_idx(t) = m;
        end
        % dummy column → coast (msmt_idx stays 0)
    end
else
    % Fallback: linear_assignment (munkres) local implementation
    [row_assign, ~] = munkres(cost_aug);
    trk_idx  = (1:num_tracks)';
    msmt_idx = zeros(num_tracks, 1);
    for ii = 1:num_tracks
        jj = row_assign(ii);
        if jj <= num_msmts && cost_aug(ii, jj) < big - 1
            msmt_idx(ii) = jj;
        end
    end
end
end


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
        %% Step 3: find an uncovered, unstarred zero
        found_zero = false;
        z_row = 0;  z_col = 0;
        for i = 1:n
            if covered_rows(i), continue; end
            for j = 1:m
                if covered_cols(j), continue; end
                if C(i,j) == 0 && ~star(i,j)
                    z_row = i;  z_col = j;
                    found_zero = true;
                    break;
                end
            end
            if found_zero, break; end
        end

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
        % Collect the full path before flipping any stars; modifying stars
        % during traversal would cause find(star(:,col)) to re-discover
        % nodes just starred, looping forever.
        curr_col = z_col;
        path_r   = z_row;
        path_c   = z_col;
        while true
            si = find(star(:, curr_col), 1);
            if isempty(si), break; end
            pj2      = find(prime(si, :), 1);
            path_r(end+1) = si;   path_c(end+1) = curr_col;  % starred zero
            path_r(end+1) = si;   path_c(end+1) = pj2;       % primed zero
            curr_col = pj2;
        end
        % Flip: star primed zeros (odd indices), unstar starred zeros (even)
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
