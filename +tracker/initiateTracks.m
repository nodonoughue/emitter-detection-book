function [new_tracks, next_track_id, new_buffer, buffer_tracks] = ...
    initiateTracks(measurements, curr_time, motion_model, msmt_model, ...
                   gate_probability, next_track_id, buffer_msmts, buffer_tracks, ...
                   target_max_velocity, target_max_acceleration, verbose)
% initiateTracks  Two-Point Initiator: seed new tentative tracks from
%                 measurements that were not associated with any existing track.
%
% [new_tracks, next_track_id, new_buffer, buffer_tracks] = ...
%     initiateTracks(measurements, curr_time, motion_model, msmt_model, ...
%                    gate_probability, next_track_id, buffer_msmts, buffer_tracks)
%
% Each unassociated measurement is first buffered as a tentative single-point
% observation.  On subsequent calls, new measurements are paired with buffered
% single-point tracks via GNN association.  Successfully paired measurements yield
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
if nargin < 11, verbose                 = false; end
verbose = false;  % verbose printing moved to runTrackerStep for divergence diagnostics

new_tracks = {};

% Track which old buffer entries should be removed after this scan.
% Mirrors Python TwoPointInitiator semantics exactly:
%   norm(s1.position) < 1.0  → zero-pos track STAYS as permanent FA sink
%   norm(s1.position) >= 1.0 → consumed after one round regardless of
%                               whether a match was found or accepted
% Initialised here so it is valid in both branches of the if below.
remove_buf_mask = false(numel(buffer_tracks), 1);

% --- 1. Try to pair new measurements with buffered single-point tracks ------
if ~isempty(buffer_tracks) && ~isempty(measurements)
    [trk_idx, msmt_idx, unmatched_msmt_idx] = ...
        tracker.associateTracks(buffer_tracks, measurements, curr_time, ...
                                motion_model, msmt_model, gate_probability, 'gnn');

    if verbose
        n_pairs = sum(msmt_idx > 0);
        fprintf('  [INIT] t=%.1fs  %d buffer tracks, %d new msmts, %d GNN pairs\n', ...
            curr_time, numel(buffer_tracks), numel(measurements), n_pairs);
    end

    for kk = 1:numel(trk_idx)
        ti = trk_idx(kk);
        mi = msmt_idx(kk);
        if mi == 0
            % No measurement matched this buffer track — remove it unconditionally.
            % One-shot buffer semantics: no permanent FA sinks.  Matches Python
            % TwoPointInitiator which always removes unmatched buffer entries so
            % zero-position tracks cannot accumulate and drain future measurements.
            remove_buf_mask(ti) = true;
            if verbose
                s1_tmp = tracker.currState(buffer_tracks{ti});
                p1 = s1_tmp.state(s1_tmp.state_space.pos_idx);
                fprintf('    buf[%d] -> no match  buf_pos=[%.1f, %.1f, %.1f] km\n', ...
                    ti, p1(1)/1e3, p1(2)/1e3, p1(3)/1e3);
            end
            continue;
        end

        % Build a two-point State with velocity estimate
        s1 = tracker.currState(buffer_tracks{ti});  % first point
        m2 = measurements{mi};                       % second measurement
        dt = m2.time - s1.time;

        if dt <= 0 || isempty(msmt_model.least_square_fun)
            if verbose
                fprintf('    buf[%d] -> msmt[%d]  REJECT dt=%.1f\n', ti, mi, dt);
            end
            continue;
        end

        % All matched buffer tracks are consumed (one-shot buffer semantics).
        remove_buf_mask(ti) = true;

        x_pos1 = s1.state(s1.state_space.pos_idx);
        if norm(x_pos1) < 1.0
            % LS failed when this buffer entry was created — skip track formation.
            if verbose
                fprintf('    buf[%d] -> msmt[%d]  REJECT x_pos1 near zero (norm=%.2f)\n', ...
                    ti, mi, norm(x_pos1));
            end
            continue;
        end

        % Estimate second-point position, seeded from first-point estimate
        [x_pos2, ~] = msmt_model.least_square_fun(m2.zeta, x_pos1);

        % Skip if the second-point LS failed (NaN, near origin, or diverged beyond 5000 km).
        % Matches Python state_from_measurement sanity check: norm > 5e6 m -> rejected.
        if any(isnan(x_pos2(:))) || norm(x_pos2) < 1.0 || norm(x_pos2) > 5e6
            if verbose
                fprintf('    buf[%d] -> msmt[%d]  REJECT x_pos2 invalid (norm=%.2f km)\n', ...
                    ti, mi, norm(x_pos2)/1e3);
                fprintf('      x_pos1=[%.1f, %.1f, %.1f] km  zeta=[%s]\n', ...
                    x_pos1(1)/1e3, x_pos1(2)/1e3, x_pos1(3)/1e3, ...
                    num2str(m2.zeta(:)', '%.1f '));
            end
            continue;
        end

        % Kinematic consistency check: the two LS-estimated positions must be
        % reachable from each other.  The bound accounts for physical motion
        % plus 1-sigma position uncertainty at each end (from the CRLB), so
        % it automatically widens in poor sensor geometry (e.g. near-coplanar
        % TDOA where σ_z can be 5–20 km) and stays tight in good geometry.
        %
        %   max_disp = 10 * (max_vel * dt
        %                    + sqrt(trace(CRLB(x_pos1)))
        %                    + sqrt(trace(CRLB(x_pos2))))
        %
        % Falls back to 10 * max_vel * dt if crlb_fun is unavailable.
        if ~isempty(target_max_velocity)
            kinematic_bound = target_max_velocity * abs(dt);
            if ~isempty(msmt_model.crlb_fun)
                crlb1 = msmt_model.crlb_fun(x_pos1);
                crlb2 = msmt_model.crlb_fun(x_pos2);
                if ~any(isnan(crlb1(:)))
                    kinematic_bound = kinematic_bound + sqrt(trace(crlb1));
                end
                if ~any(isnan(crlb2(:)))
                    kinematic_bound = kinematic_bound + sqrt(trace(crlb2));
                end
            end
            max_disp = 10.0 * kinematic_bound;
            if norm(x_pos2 - x_pos1) > max_disp
                if verbose
                    fprintf('    buf[%d] -> msmt[%d]  REJECT kinematic inconsistency (disp=%.1f km > %.1f km)\n', ...
                        ti, mi, norm(x_pos2-x_pos1)/1e3, max_disp/1e3);
                end
                continue;
            end
        end

        % Velocity estimate: finite difference between LS-estimated positions.
        % Clip to target_max_velocity so the first prediction step does not
        % move the track to an unphysical position before constrainMotion fires.
        vel_est = (x_pos2 - x_pos1) / dt;
        if ~isempty(target_max_velocity) && norm(vel_est) > target_max_velocity
            vel_est = vel_est * (target_max_velocity / norm(vel_est));
        end

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

        if verbose
            fprintf('    buf[%d] -> msmt[%d]  ACCEPT  new track_id=%d\n', ...
                ti, mi, new_trk.track_id);
            fprintf('      x_pos1=[%.1f, %.1f, %.1f] km\n', ...
                x_pos1(1)/1e3, x_pos1(2)/1e3, x_pos1(3)/1e3);
            fprintf('      x_pos2=[%.1f, %.1f, %.1f] km  dt=%.1fs\n', ...
                x_pos2(1)/1e3, x_pos2(2)/1e3, x_pos2(3)/1e3, dt);
            fprintf('      vel_est=[%.1f, %.1f, %.1f] m/s (speed=%.1f)\n', ...
                vel_est(1), vel_est(2), vel_est(3), norm(vel_est));
        end

    end

    % Only truly unmatched measurements (not assigned by GNN at all) become new
    % buffer entries.  Measurements paired with a rejected buffer track are consumed
    % without recycling — matching Python TwoPointInitiator behaviour where rejected
    % buffer tracks stay as permanent FA-measurement sinks.
    unmatched_new = measurements(unmatched_msmt_idx);
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
% Retain old buffer entries that were NOT consumed: zero-position sinks
% (norm < 1 m) stay indefinitely; non-zero tracks are removed after one
% pairing attempt whether or not it produced a confirmed track.
new_buffer    = buffer_msmts(~remove_buf_mask);
buffer_tracks = buffer_tracks(~remove_buf_mask);

for jj = 1:numel(unmatched_new)
    m = unmatched_new{jj};

    ss    = motion_model.state_space;
    n_pos = numel(ss.pos_idx);
    x_buf = zeros(ss.num_states, 1);
    P_buf = 1e6 * eye(ss.num_states);  % fallback: wide uncertainty

    if ~isempty(msmt_model.least_square_fun)
        try
            [x_pos_est, ~] = msmt_model.least_square_fun(m.zeta, zeros(n_pos, 1));
            if any(isnan(x_pos_est(:)))
                error('tracker:initiateTracks:lsNaN', 'LS returned NaN; falling back to x=0');
            end

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

    if verbose
        fprintf('  [BUF+] jj=%d  x_buf=[%.1f, %.1f, %.1f] km (norm=%.1f km)  zeta=[%s]\n', ...
            jj, x_buf(ss.pos_idx(1))/1e3, x_buf(ss.pos_idx(2))/1e3, x_buf(ss.pos_idx(3))/1e3, ...
            norm(x_buf(ss.pos_idx))/1e3, num2str(m.zeta(:)', '%.1f '));
        P_pos_diag = diag(P_buf(ss.pos_idx, ss.pos_idx));
        fprintf('         P_pos_diag=[%.2e, %.2e, %.2e] m^2\n', ...
            P_pos_diag(1), P_pos_diag(2), P_pos_diag(3));
    end
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
