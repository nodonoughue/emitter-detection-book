function [new_tracks, next_track_id, new_buffer, buffer_tracks] = ...
    initiateTracks(measurements, curr_time, motion_model, msmt_model, ...
                   gate_probability, next_track_id, buffer_msmts, buffer_tracks, ...
                   target_max_velocity, target_max_acceleration)
% initiateTracks  Two-Point Initiator: seed new tentative tracks from
%                 measurements that were not associated with any existing track.
%
% [new_tracks, next_track_id, new_buffer, buffer_tracks] = ...
%     initiateTracks(measurements, curr_time, motion_model, msmt_model, ...
%                    gate_probability, next_track_id, buffer_msmts, buffer_tracks)
%
% Each unassociated measurement is first buffered as a tentative single-point
% observation.  On subsequent calls, new measurements are paired with buffered
% single-point tracks via NN association.  Successfully paired measurements yield
% a new tentative two-point Track (with a direct velocity estimate).  Unmatched
% new measurements replace the buffer for the next scan.
%
% INPUTS
%   measurements    Cell array of unassociated Measurement structs
%   curr_time       Timestamp of the current scan [s]
%   motion_model    Motion model struct from makeMotionModel
%   msmt_model      Measurement model struct from makeMeasurementModel
%                   (must have a non-empty .least_square_fun field)
%   gate_probability         Gate probability for buffer–measurement pairing
%   next_track_id            Integer: ID to assign to the next new tentative track
%   buffer_msmts             Cell array of Measurement structs buffered from the
%                            previous scan (pass {} on the first call)
%   buffer_tracks            Cell array of single-point Track structs corresponding
%                            to buffer_msmts (pass {} on the first call)
%   target_max_velocity      Optional max target speed [m/s].  Caps the velocity
%                            covariance block and stamps max_velocity on each new
%                            track so runTrackerStep can clip the EKF state.
%                            Pass [] to disable (default: [])
%   target_max_acceleration  Optional max target acceleration [m/s²].  Caps the
%                            acceleration covariance block and stamps max_acceleration
%                            on each new track.  Pass [] to disable (default: [])
%
% OUTPUTS
%   new_tracks      Cell array of newly created tentative two-point Track structs
%   next_track_id   Updated track ID counter
%   new_buffer      Cell array of Measurement structs to carry into next call
%   buffer_tracks   Cell array of Track structs corresponding to new_buffer
%
% Nicholas O'Donoughue
% June 2025

if nargin < 9,  target_max_velocity     = []; end
if nargin < 10, target_max_acceleration = []; end

new_tracks = {};

% --- 1. Try to pair new measurements with buffered single-point tracks ------
if ~isempty(buffer_tracks) && ~isempty(measurements)
    [trk_idx, msmt_idx, unmatched_msmt_idx] = ...
        tracker.associateTracks(buffer_tracks, measurements, curr_time, ...
                                motion_model, msmt_model, gate_probability, 'nn');

    matched_new_msmt_mask = false(numel(measurements), 1);

    for kk = 1:numel(trk_idx)
        ti = trk_idx(kk);
        mi = msmt_idx(kk);
        if mi == 0
            continue;  % no match for this buffered track
        end

        % Build a two-point State with velocity estimate
        s1 = tracker.currState(buffer_tracks{ti});  % first point
        m2 = measurements{mi};                       % second measurement
        dt = m2.time - s1.time;

        if dt <= 0 || isempty(msmt_model.least_square_fun)
            continue;
        end

        % Use the position already stored in the buffer track.  Buffer tracks
        % are now initialised via LS (with a try-catch + 500 km sanity check),
        % so s1.state(pos_idx) is the best available first-point estimate.
        % Re-running LS from zeros here can diverge (targets are 50-125 km
        % from the origin; the divergence detector fires before convergence),
        % producing vel_est = x_pos2/dt — thousands of km/s — and causing
        % the track to predict wildly off-course after just one step.
        % This matches Python TwoPointInitiator, which uses s1.position directly.
        x_pos1 = s1.state(s1.state_space.pos_idx);
        if norm(x_pos1) < 1.0
            continue;   % buffer LS failed (fell back to zeros); skip this pair
        end

        % Estimate second-point position, seeded from first-point estimate
        [x_pos2, ~] = msmt_model.least_square_fun(m2.zeta, x_pos1);

        % Velocity estimate: finite difference between LS-estimated positions
        vel_est = (x_pos2 - x_pos1) / dt;

        % Build initial state vector
        ss = motion_model.state_space;
        x_init = zeros(ss.num_states, 1);
        x_init(ss.pos_idx) = x_pos2;
        if ss.has_vel
            x_init(ss.vel_idx) = vel_est;
        end

        % Initial covariance: use CRLB from measurement Jacobian (10x for conservatism),
        % scaled by dt^2 / dt^4 for velocity / acceleration blocks.
        % This matches Python TwoPointInitiator._build_initial_covariance.
        P_init = build_initial_covariance(ss, x_pos2, vel_est, m2.time, dt, msmt_model, ...
                                          target_max_velocity, target_max_acceleration);

        s_init = tracker.makeState(ss, m2.time, x_init, P_init);
        new_trk = tracker.makeTrack(s_init, next_track_id);
        new_trk.max_velocity    = target_max_velocity;
        new_trk.max_acceleration = target_max_acceleration;
        next_track_id = next_track_id + 1;
        new_tracks{end+1} = new_trk;  %#ok<AGROW>

        matched_new_msmt_mask(mi) = true;
    end

    % Measurements that were not used to form two-point tracks become the new buffer
    unmatched_new = measurements(~matched_new_msmt_mask);
else
    % No buffered tracks yet — all new measurements go to buffer
    unmatched_new = measurements;
end

% --- 2. Build single-point buffer tracks from unmatched new measurements ---
% We estimate the position for each buffer track via LS so that the gate is
% centred on the estimated source location.  Without this, every buffer track
% has position = zeros and the NN pairing at the next scan is random — a
% true-target buffer track can be paired with a false-alarm new measurement,
% producing a completely wrong velocity even after the two-point fix.
% LS is run with a try-catch; infeasible false-alarm measurements cause the
% solver to diverge, and we fall back to x=0 / P=1e6*I for those.
new_buffer    = {};
buffer_tracks = {};

for jj = 1:numel(unmatched_new)
    m = unmatched_new{jj};

    ss    = motion_model.state_space;
    n_pos = numel(ss.pos_idx);
    x_buf = zeros(ss.num_states, 1);
    P_buf = 1e6 * eye(ss.num_states);  % fallback: wide uncertainty

    if ~isempty(msmt_model.least_square_fun)
        try
            [x_pos_est, ~] = msmt_model.least_square_fun(m.zeta, zeros(n_pos, 1));

            % If 3D and z < 0, retry with the z sign flipped to escape the
            % mirror-image local minimum from near-coplanar sensor arrays.
            if n_pos >= 3 && x_pos_est(3) < 0
                x_retry = zeros(n_pos, 1);
                x_retry(3) = -x_pos_est(3);
                [x_pos_est_retry, ~] = msmt_model.least_square_fun(m.zeta, x_retry);
                if x_pos_est_retry(3) >= 0
                    x_pos_est = x_pos_est_retry;
                end
            end

            % Accept the LS result only if it looks physically reasonable
            % (within 500 km of the origin; avoids using a diverged iterate)
            if norm(x_pos_est) < 5e5
                x_buf(ss.pos_idx) = x_pos_est;

                % Position covariance: use the dedicated CRLB function when
                % available (numerically robust pinv; handles near-singular
                % geometries from false-alarm LS solutions).  Fall back to
                % manual H'R^{-1}H FIM only when crlb_fun is not provided.
                if ~isempty(msmt_model.crlb_fun)
                    crlb_pos = msmt_model.crlb_fun(x_pos_est);
                    if ~any(isnan(crlb_pos(:)))
                        P_buf(ss.pos_idx, ss.pos_idx) = crlb_pos;
                    end
                else
                    s_tmp = tracker.makeState(ss, m.time, x_buf, eye(ss.num_states));
                    H_raw = msmt_model.h_fun(s_tmp);
                    if ndims(H_raw) == 3, H_raw = H_raw(:,:,1); end
                    H_pos = H_raw(:, ss.pos_idx);
                    FIM   = H_pos' * (msmt_model.R \ H_pos);
                    if rcond(FIM) >= 1e-12
                        P_buf(ss.pos_idx, ss.pos_idx) = 10.0 * inv(FIM);
                    end
                    % Velocity/accel blocks: no time-difference available yet;
                    % keep at 1e6 (large uncertainty — corrected at two-point stage)
                end
            end
        catch
            % LS failed; keep x=0, P=1e6*I
        end
    end

    % Apply velocity/acceleration covariance caps to buffer track P
    if ~isempty(target_max_velocity) && ss.has_vel && ~isempty(ss.vel_idx)
        vel_block = P_buf(ss.vel_idx, ss.vel_idx);
        max_diag  = max(diag(vel_block));
        if max_diag > target_max_velocity^2
            P_buf(ss.vel_idx, ss.vel_idx) = vel_block * (target_max_velocity^2 / max_diag);
        end
    end
    if ~isempty(target_max_acceleration) && isfield(ss, 'accel_idx') && ~isempty(ss.accel_idx)
        acc_block = P_buf(ss.accel_idx, ss.accel_idx);
        max_diag  = max(diag(acc_block));
        if max_diag > target_max_acceleration^2
            P_buf(ss.accel_idx, ss.accel_idx) = acc_block * (target_max_acceleration^2 / max_diag);
        end
    end

    s_buf  = tracker.makeState(ss, m.time, x_buf, P_buf);
    t_buf  = tracker.makeTrack(s_buf, -1);  % id=-1: tentative buffer track
    t_buf.max_velocity    = target_max_velocity;
    t_buf.max_acceleration = target_max_acceleration;
    new_buffer{end+1}    = m;     %#ok<AGROW>
    buffer_tracks{end+1} = t_buf; %#ok<AGROW>
end


%% ---- Local helper -----------------------------------------------------------

function P = build_initial_covariance(ss, x_pos, vel_est, t, dt, msmt_model, ...
                                       target_max_velocity, target_max_acceleration)
% Compute the initial state covariance for a two-point tentative track using
% a 10x Fisher-information CRLB, matching Python TwoPointInitiator logic.
%
% Position block  : 10 * crlb_pos  (via crlb_fun or manual FIM)
% Velocity block  : crlb_pos / dt^2
% Accel block     : crlb_pos / dt^4  (CA/CJ only)
%
% Falls back to 1e6*I if the CRLB is unavailable or ill-conditioned.

if nargin < 7, target_max_velocity     = []; end
if nargin < 8, target_max_acceleration = []; end

pos_covar_multiplier = 10.0;

% Prefer the dedicated CRLB function (uses pinv; robust to near-singular
% geometries).  Fall back to the manual H'R^{-1}H inversion when absent.
if ~isempty(msmt_model.crlb_fun)
    crlb_raw = msmt_model.crlb_fun(x_pos);   % (num_dims x num_dims)
    if any(isnan(crlb_raw(:)))
        P = 1e6 * eye(ss.num_states);
        return;
    end
    crlb_pos = pos_covar_multiplier * crlb_raw;
else
    % Build a minimal state struct at the estimated position
    x_temp = zeros(ss.num_states, 1);
    x_temp(ss.pos_idx) = x_pos;
    if ss.has_vel && ~isempty(ss.vel_idx)
        x_temp(ss.vel_idx) = vel_est;
    end
    s_temp = tracker.makeState(ss, t, x_temp, eye(ss.num_states));

    H = msmt_model.h_fun(s_temp);
    if ndims(H) == 3
        H = H(:, :, 1)';
    end
    H_pos = H(:, ss.pos_idx);
    FIM   = H_pos' * (msmt_model.R \ H_pos);

    if rcond(FIM) < 1e-12
        P = 1e6 * eye(ss.num_states);
        return;
    end
    crlb_pos = pos_covar_multiplier * inv(FIM);
end

crlb_vel   = crlb_pos / (dt^2);
crlb_accel = crlb_vel / (dt^2);

% Apply velocity/acceleration covariance caps if requested
if ~isempty(target_max_velocity) && ss.has_vel && ~isempty(ss.vel_idx)
    max_diag = max(diag(crlb_vel));
    if max_diag > target_max_velocity^2
        crlb_vel = crlb_vel * (target_max_velocity^2 / max_diag);
    end
end
if ~isempty(target_max_acceleration) && isfield(ss, 'accel_idx') && ~isempty(ss.accel_idx)
    max_diag = max(diag(crlb_accel));
    if max_diag > target_max_acceleration^2
        crlb_accel = crlb_accel * (target_max_acceleration^2 / max_diag);
    end
end

P = zeros(ss.num_states);
P(ss.pos_idx, ss.pos_idx) = crlb_pos;
if ss.has_vel && ~isempty(ss.vel_idx)
    P(ss.vel_idx, ss.vel_idx) = crlb_vel;
end
if isfield(ss, 'accel_idx') && ~isempty(ss.accel_idx)
    P(ss.accel_idx, ss.accel_idx) = crlb_accel;
end
