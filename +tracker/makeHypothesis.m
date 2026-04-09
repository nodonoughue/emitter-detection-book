function h = makeHypothesis(track, msmt, motion_model, curr_time, gate_probability, num_msmt_dims)
% makeHypothesis  Create a hypothesis struct associating a track with a measurement.
%
% Normal form — associate a measurement with a track:
%   h = makeHypothesis(track, msmt, motion_model)
%
% Missed-detection form — coast a track with no measurement:
%   h = makeHypothesis(track, [], motion_model, curr_time)
%   h = makeHypothesis(track, [], motion_model, curr_time, gate_probability, num_msmt_dims)
%
% INPUTS (normal form)
%   track        Track struct from makeTrack
%   msmt         Measurement struct from makeMeasurement (must have .msmt_model set)
%   motion_model Motion model struct from makeMotionModel
%
% INPUTS (missed-detection form)
%   track            Track struct from makeTrack
%   msmt             [] (empty signals a missed detection)
%   motion_model     Motion model struct from makeMotionModel
%   curr_time        Timestamp to coast the track to [s]
%   gate_probability Chi-square gate probability (e.g., 0.99); combined with
%                    num_msmt_dims to compute distance = chi2inv(gate_probability, num_msmt_dims).
%                    Optional: if omitted, distance defaults to 0.0 (coasting always wins,
%                    appropriate when there are no competing measurements in the scan).
%   num_msmt_dims    Measurement vector dimension; required alongside gate_probability to
%                    compute the gate threshold. Omit when no measurements are available.
%
% OUTPUTS
%   h  Hypothesis struct with fields:
%        track         — copy of the input track struct
%        msmt          — copy of the input measurement struct (or [])
%        motion_model  — copy of the motion model struct
%        s_pred        — predicted State struct (track coasted to msmt.time or curr_time)
%        z_pred        — predicted measurement vector (num_msmt x 1), or [] for missed
%        H             — measurement Jacobian (num_msmt x num_states), or [] for missed
%        innov         — innovation: msmt.zeta - z_pred (num_msmt x 1), or [] for missed
%        S             — innovation covariance H*P*H' + R (num_msmt x num_msmt), or [] for missed
%        distance      — Mahalanobis distance squared: innov' * S^{-1} * innov
%        likelihood    — multivariate Gaussian likelihood N(innov; 0, S)
%        log_likelihood — log of likelihood
%        is_valid      — true; set to false by apply_gate to mark out-of-gate hypotheses
%        is_missed     — true for missed-detection hypotheses, false for normal ones
%
% Gating
%   To apply a chi-square acceptance gate after construction:
%       gate_size = chi2inv(gate_probability, numel(h.innov));
%       if h.distance > gate_size
%           h.is_valid    = false;
%           h.distance    = Inf;
%           h.likelihood  = 0;
%           h.log_likelihood = -Inf;
%       end
%
% Missed-detection hypotheses
%   The missed-detection form uses the same struct layout, making it easy to
%   mix normal and missed-detection hypotheses in a cost matrix or likelihood
%   table without special-casing.  The is_missed flag identifies them.
%   Fields z_pred, H, innov, and S are all [] for missed detections.
%
% Nicholas O'Donoughue
% June 2025

is_missed = isempty(msmt);

if is_missed
    %% --- Missed-detection hypothesis -----------------------------------
    if nargin < 4 || isempty(curr_time)
        error('tracker:makeHypothesis:missingTime', ...
              'curr_time is required for missed-detection hypotheses (msmt = []).');
    end

    if nargin >= 6 && ~isempty(gate_probability) && ~isempty(num_msmt_dims)
        % Self-infer: gate threshold chi2inv(gate_probability, num_msmt_dims)
        distance = chi2inv(gate_probability, num_msmt_dims);
    else
        % No competing measurements in this scan; coasting always wins
        distance = 0.0;
    end

    s      = tracker.currState(track);
    s_pred = tracker.predictState(s, curr_time, motion_model);

    h = struct( ...
        'track',         track,      ...
        'msmt',          [],         ...
        'motion_model',  motion_model, ...
        's_pred',        s_pred,     ...
        'z_pred',        [],         ...
        'H',             [],         ...
        'innov',         [],         ...
        'S',             [],         ...
        'distance',      distance,   ...
        'likelihood',    distance,   ...   % re-use distance as weight (PDA / NN convention)
        'log_likelihood', log(max(distance, realmin)), ...
        'is_valid',      true,       ...
        'is_missed',     true);

else
    %% --- Normal hypothesis (track + measurement) -----------------------
    if ~isfield(msmt, 'msmt_model') || isempty(msmt.msmt_model)
        error('tracker:makeHypothesis:missingMsmtModel', ...
              'msmt.msmt_model must be set. Create measurements with makeMeasurement(msmt_model, state, time).');
    end

    msmt_model = msmt.msmt_model;

    % Predict track to measurement time
    s      = tracker.currState(track);
    s_pred = tracker.predictState(s, msmt.time, motion_model);

    % Predicted measurement and Jacobian
    z_pred = msmt_model.z_fun(s_pred);
    H      = msmt_model.h_fun(s_pred);

    % Squeeze 3-D Jacobian (num_states x num_msmt x 1) -> (num_msmt x num_states)
    if ndims(H) == 3
        H = H(:, :, 1)';
    end

    % Innovation and innovation covariance
    innov = msmt.zeta(:) - z_pred(:);
    S     = H * s_pred.covar * H' + msmt_model.R;

    % Mahalanobis distance squared: innov' * S^{-1} * innov
    dist = innov' / S * innov;

    % Multivariate Gaussian likelihood: N(innov; 0, S)
    n       = numel(innov);
    log_lik = -0.5 * (n * log(2*pi) + log(max(det(S), realmin)) + dist);
    lik     = exp(log_lik);

    h = struct( ...
        'track',         track,      ...
        'msmt',          msmt,       ...
        'motion_model',  motion_model, ...
        's_pred',        s_pred,     ...
        'z_pred',        z_pred(:),  ...
        'H',             H,          ...
        'innov',         innov,      ...
        'S',             S,          ...
        'distance',      dist,       ...
        'likelihood',    lik,        ...
        'log_likelihood', log_lik,   ...
        'is_valid',      true,       ...
        'is_missed',     false);
end
