function [new_tracks, next_track_id, new_buf1_msmts, new_buf1_tracks, ...
          new_buf2_msmts, new_buf2_tracks] = ...
    initiateTracks(measurements, curr_time, motion_model, msmt_model, ...
                   gate_probability, next_track_id, buf1_msmts, buf1_tracks, ...
                   target_max_velocity, target_max_acceleration, verbose, ...
                   init_type, buf2_msmts, buf2_tracks)
% initiateTracks  Initiate new tentative tracks from unassociated measurements.
%
% Three initiation strategies are available:
%
%   'one_point'    SinglePointMeasurementInitiator — one tentative track per
%                  measurement, position from LS only (velocity initialised
%                  to zero with large uncertainty).  No buffering.
%
%   'two_point'    TwoPointInitiator (default) — buffers single-point
%                  observations and pairs them across two consecutive scans
%                  to obtain a velocity estimate via finite difference.
%
%   'three_point'  ThreePointInitiator — buffers single-point (stage 1) and
%                  two-point (stage 2) observations across three consecutive
%                  scans to obtain velocity AND acceleration estimates via
%                  central finite differences.
%
% [new_tracks, next_track_id, new_buf1_msmts, new_buf1_tracks,           ...
%  new_buf2_msmts, new_buf2_tracks] =                                     ...
%     initiateTracks(measurements, curr_time, motion_model, msmt_model,   ...
%                    gate_probability, next_track_id,                      ...
%                    buf1_msmts, buf1_tracks,                              ...
%                    target_max_velocity, target_max_acceleration, verbose, ...
%                    init_type, buf2_msmts, buf2_tracks)
%
% INPUTS
%   measurements    Cell array of unassociated Measurement structs
%   curr_time       Timestamp of the current scan [s]
%   motion_model    Motion model struct from makeMotionModel
%   msmt_model      Measurement model struct from makeMeasurementModel
%                   (must have a non-empty .least_square_fun field)
%   gate_probability         Gate probability for buffer-measurement pairing
%   next_track_id            Integer: ID to assign to the next new track
%   buf1_msmts               Stage-1 measurement buffer (pass {} on first call)
%   buf1_tracks              Stage-1 track buffer (pass {} on first call)
%   target_max_velocity      Optional max target speed [m/s] (default: [])
%   target_max_acceleration  Optional max target acceleration [m/s²] (default: [])
%   verbose                  Print diagnostic information (default: false)
%   init_type                'one_point', 'two_point' (default), or 'three_point'
%   buf2_msmts               Stage-2 measurement buffer; only used by
%                            'three_point' (pass {} otherwise)
%   buf2_tracks              Stage-2 track buffer; only used by 'three_point'
%
% OUTPUTS
%   new_tracks       Cell array of newly created tentative Track structs
%   next_track_id    Updated track ID counter
%   new_buf1_msmts   Updated stage-1 measurement buffer
%   new_buf1_tracks  Updated stage-1 track buffer
%   new_buf2_msmts   Updated stage-2 measurement buffer (empty for one/two_point)
%   new_buf2_tracks  Updated stage-2 track buffer (empty for one/two_point)
%
% Nicholas O'Donoughue
% June 2025

if nargin < 9,  target_max_velocity     = []; end
if nargin < 10, target_max_acceleration = []; end
if nargin < 11, verbose                 = false; end
if nargin < 12 || isempty(init_type),   init_type  = 'two_point'; end
if nargin < 13, buf2_msmts  = {}; end
if nargin < 14, buf2_tracks = {}; end

% Passthrough defaults for unused stage-2 buffers
new_buf2_msmts  = {};
new_buf2_tracks = {};

switch lower(init_type)
    case {'one_point', 'single_point'}
        [new_tracks, next_track_id] = ...
            initiate_one_point(measurements, motion_model, msmt_model, ...
                               target_max_velocity, target_max_acceleration, ...
                               next_track_id, verbose);
        new_buf1_msmts  = buf1_msmts;
        new_buf1_tracks = buf1_tracks;

    case 'two_point'
        [new_tracks, next_track_id, new_buf1_msmts, new_buf1_tracks] = ...
            initiate_two_point(measurements, curr_time, motion_model, msmt_model, ...
                               gate_probability, next_track_id, ...
                               buf1_msmts, buf1_tracks, ...
                               target_max_velocity, target_max_acceleration, verbose);

    case 'three_point'
        [new_tracks, next_track_id, new_buf1_msmts, new_buf1_tracks, ...
         new_buf2_msmts, new_buf2_tracks] = ...
            initiate_three_point(measurements, curr_time, motion_model, msmt_model, ...
                                 gate_probability, next_track_id, ...
                                 buf1_msmts, buf1_tracks, ...
                                 buf2_msmts, buf2_tracks, ...
                                 target_max_velocity, target_max_acceleration, verbose);

    otherwise
        error('tracker:initiateTracks:unknownType', ...
              'Unknown init_type "%s". Use ''one_point'', ''two_point'', or ''three_point''.', ...
              init_type);
end

end  % initiateTracks


%% ==========================================================================
function [new_tracks, next_track_id] = ...
    initiate_one_point(measurements, motion_model, msmt_model, ...
                       target_max_velocity, target_max_acceleration, ...
                       next_track_id, verbose)
% SinglePointMeasurementInitiator: one tentative track per measurement.
% No buffering; velocity is initialised to zero with large uncertainty.

new_tracks = {};

for jj = 1:numel(measurements)
    [t_new, ok] = build_single_point_track(measurements{jj}, motion_model, msmt_model, ...
                                           target_max_velocity, target_max_acceleration, ...
                                           next_track_id);
    if ~ok
        if verbose
            fprintf('  [1PT skip] jj=%d  LS failed — measurement discarded\n', jj);
        end
        continue;
    end
    next_track_id = next_track_id + 1;
    new_tracks{end+1} = t_new;  %#ok<AGROW>
    if verbose
        ss  = t_new.states{1}.state_space;
        pos = t_new.states{1}.state(ss.pos_idx);
        fprintf('  [1PT new ] jj=%d  track_id=%d  pos=[%.1f %.1f %.1f] km\n', ...
            jj, t_new.track_id, pos(1)/1e3, pos(2)/1e3, pos(3)/1e3);
    end
end

end  % initiate_one_point


%% ==========================================================================
function [new_tracks, next_track_id, new_buf1_msmts, new_buf1_tracks] = ...
    initiate_two_point(measurements, curr_time, motion_model, msmt_model, ...
                       gate_probability, next_track_id, buf1_msmts, buf1_tracks, ...
                       target_max_velocity, target_max_acceleration, verbose)
% TwoPointInitiator: pairs measurements across two consecutive scans for a
% direct velocity estimate via finite difference.

new_tracks = {};
remove_buf_mask = false(numel(buf1_tracks), 1);

%% Step 1: pair new measurements with buffered single-point tracks
if ~isempty(buf1_tracks) && ~isempty(measurements)
    [trk_idx, msmt_idx, unmatched_msmt_idx] = ...
        tracker.associateTracks(buf1_tracks, measurements, curr_time, ...
                                motion_model, msmt_model, gate_probability, 'gnn');

    if verbose
        fprintf('  [INIT 2PT] t=%.1fs  %d buf, %d msmts, %d pairs\n', ...
            curr_time, numel(buf1_tracks), numel(measurements), sum(msmt_idx > 0));
    end

    for kk = 1:numel(trk_idx)
        ti = trk_idx(kk);
        mi = msmt_idx(kk);

        if mi == 0
            % Unmatched buffer track — remove unconditionally (one-shot)
            remove_buf_mask(ti) = true;
            continue;
        end

        s1 = tracker.currState(buf1_tracks{ti});
        m2 = measurements{mi};
        dt = m2.time - s1.time;

        if dt <= 0 || isempty(msmt_model.least_square_fun)
            continue;
        end

        % All matched buffer tracks are consumed (one-shot)
        remove_buf_mask(ti) = true;

        x_pos1 = s1.state(s1.state_space.pos_idx);
        if norm(x_pos1) < 1.0
            continue;  % LS failed when buffer track was created
        end

        [x_pos2, ~] = msmt_model.least_square_fun(m2.zeta, x_pos1);
        if any(isnan(x_pos2(:))) || norm(x_pos2) < 1.0 || norm(x_pos2) > 5e6
            continue;
        end

        % Kinematic consistency check
        if ~isempty(target_max_velocity)
            kinematic_bound = target_max_velocity * abs(dt);
            if ~isempty(msmt_model.crlb_fun)
                crlb1 = msmt_model.crlb_fun(x_pos1);
                crlb2 = msmt_model.crlb_fun(x_pos2);
                if ~any(isnan(crlb1(:))), kinematic_bound = kinematic_bound + sqrt(trace(crlb1)); end
                if ~any(isnan(crlb2(:))), kinematic_bound = kinematic_bound + sqrt(trace(crlb2)); end
            end
            if norm(x_pos2 - x_pos1) > 10.0 * kinematic_bound
                continue;
            end
        end

        vel_est = (x_pos2 - x_pos1) / dt;
        if ~isempty(target_max_velocity) && norm(vel_est) > target_max_velocity
            vel_est = vel_est * (target_max_velocity / norm(vel_est));
        end

        ss = motion_model.state_space;
        x_init = zeros(ss.num_states, 1);
        x_init(ss.pos_idx) = x_pos2;
        if ss.has_vel
            x_init(ss.vel_idx) = vel_est;
        end

        P_init = build_initial_covariance_two_point(ss, x_pos2, vel_est, m2.time, dt, ...
                                                    msmt_model, target_max_velocity, ...
                                                    target_max_acceleration);
        s_init  = tracker.makeState(ss, m2.time, x_init, P_init);
        new_trk = tracker.makeTrack(s_init, next_track_id);
        new_trk.max_velocity     = target_max_velocity;
        new_trk.max_acceleration = target_max_acceleration;
        next_track_id = next_track_id + 1;
        new_tracks{end+1} = new_trk;  %#ok<AGROW>

        if verbose
            fprintf('  [2PT ACCEPT] buf[%d]+msmt[%d] -> track_id=%d  pos=[%.1f %.1f %.1f] km  spd=%.1f m/s\n', ...
                ti, mi, new_trk.track_id, ...
                x_pos2(1)/1e3, x_pos2(2)/1e3, x_pos2(3)/1e3, norm(vel_est));
        end
    end

    unmatched_new = measurements(unmatched_msmt_idx);
else
    unmatched_new = measurements;
end

%% Step 2: build stage-1 buffer from unmatched new measurements
new_buf1_msmts  = buf1_msmts(~remove_buf_mask);
new_buf1_tracks = buf1_tracks(~remove_buf_mask);

for jj = 1:numel(unmatched_new)
    [t_buf, ok] = build_single_point_track(unmatched_new{jj}, motion_model, msmt_model, ...
                                           target_max_velocity, target_max_acceleration, -1);
    if ~ok
        continue;  % LS failed; discard
    end
    new_buf1_msmts{end+1}  = unmatched_new{jj};  %#ok<AGROW>
    new_buf1_tracks{end+1} = t_buf;               %#ok<AGROW>
end

end  % initiate_two_point


%% ==========================================================================
function [new_tracks, next_track_id, new_buf1_msmts, new_buf1_tracks, ...
          new_buf2_msmts, new_buf2_tracks] = ...
    initiate_three_point(measurements, curr_time, motion_model, msmt_model, ...
                         gate_probability, next_track_id, ...
                         buf1_msmts, buf1_tracks, ...
                         buf2_msmts, buf2_tracks, ...
                         target_max_velocity, target_max_acceleration, verbose)
% ThreePointInitiator: buffers measurements across three consecutive scans.
% Confirmed tracks carry position, velocity, and (when the state space
% includes acceleration) acceleration estimates via central finite differences.

new_tracks    = {};
new_buf2_msmts  = {};
new_buf2_tracks = {};

%% Step 1: advance stage-2 tracks (t1, t2) to confirmed tracks using s3
unmatched_from_2 = measurements;

if ~isempty(buf2_tracks) && ~isempty(measurements)
    [trk_idx2, msmt_idx2, unmatched_msmt_idx2] = ...
        tracker.associateTracks(buf2_tracks, measurements, curr_time, ...
                                motion_model, msmt_model, gate_probability, 'gnn');

    if verbose
        fprintf('  [INIT 3PT] t=%.1fs  stage-2 %d buf, %d msmts, %d pairs\n', ...
            curr_time, numel(buf2_tracks), numel(measurements), sum(msmt_idx2 > 0));
    end

    for kk = 1:numel(trk_idx2)
        ti2 = trk_idx2(kk);
        mi2 = msmt_idx2(kk);

        if mi2 == 0
            continue;  % missed detection at stage 2 — retire without promotion
        end

        trk2   = buf2_tracks{ti2};
        s1     = trk2.states{1};              % first (stage-1) state
        s2     = tracker.currState(trk2);     % second (stage-2) state
        m3     = measurements{mi2};

        x_pos1 = s1.state(s1.state_space.pos_idx);
        x_pos2 = s2.state(s2.state_space.pos_idx);

        % Third-point LS, seeded from stage-2 position
        [x_pos3, ~] = msmt_model.least_square_fun(m3.zeta, x_pos2);
        if any(isnan(x_pos3(:))) || norm(x_pos3) < 1.0 || norm(x_pos3) > 5e6
            continue;
        end

        dt1 = s2.time - s1.time;
        dt2 = m3.time - s2.time;

        % Require consistent time steps (within 10%)
        if dt1 <= 0 || dt2 <= 0 || abs(dt1 - dt2) / max(dt1, dt2) > 0.1
            continue;
        end
        dt = (dt1 + dt2) / 2.0;

        % Finite-difference estimates
        vel_est = (x_pos3 - x_pos1) / (2.0 * dt);   % central difference
        acc_est = (x_pos3 - 2*x_pos2 + x_pos1) / (dt^2);  % second difference

        if ~isempty(target_max_velocity) && norm(vel_est) > target_max_velocity
            vel_est = vel_est * (target_max_velocity / norm(vel_est));
        end

        ss = motion_model.state_space;
        x_init = zeros(ss.num_states, 1);
        x_init(ss.pos_idx) = x_pos3;
        if ss.has_vel
            x_init(ss.vel_idx) = vel_est;
        end
        if isfield(ss, 'accel_idx') && ~isempty(ss.accel_idx)
            x_init(ss.accel_idx) = acc_est;
        end

        P_init  = build_initial_covariance_three_point(ss, x_pos1, x_pos2, x_pos3, dt, ...
                                                       msmt_model, target_max_velocity, ...
                                                       target_max_acceleration);
        s_init  = tracker.makeState(ss, m3.time, x_init, P_init);
        new_trk = tracker.makeTrack(s_init, trk2.track_id);  % inherit track ID
        new_trk.max_velocity     = target_max_velocity;
        new_trk.max_acceleration = target_max_acceleration;
        new_tracks{end+1} = new_trk;  %#ok<AGROW>

        if verbose
            fprintf('  [3PT ACCEPT] trk_id=%d  pos=[%.1f %.1f %.1f] km  spd=%.1f m/s\n', ...
                trk2.track_id, x_pos3(1)/1e3, x_pos3(2)/1e3, x_pos3(3)/1e3, norm(vel_est));
        end
    end

    unmatched_from_2 = measurements(unmatched_msmt_idx2);
end
% Stage-2 tracks are ALL retired one-shot (new_buf2 starts empty above)

%% Step 2: advance stage-1 tracks to stage-2 using the remaining measurements
unmatched_from_1 = unmatched_from_2;

if ~isempty(buf1_tracks) && ~isempty(unmatched_from_2)
    [trk_idx1, msmt_idx1, unmatched_msmt_idx1] = ...
        tracker.associateTracks(buf1_tracks, unmatched_from_2, curr_time, ...
                                motion_model, msmt_model, gate_probability, 'gnn');

    if verbose
        fprintf('  [INIT 3PT] t=%.1fs  stage-1 %d buf, %d msmts, %d pairs\n', ...
            curr_time, numel(buf1_tracks), numel(unmatched_from_2), sum(msmt_idx1 > 0));
    end

    for kk = 1:numel(trk_idx1)
        ti1 = trk_idx1(kk);
        mi1 = msmt_idx1(kk);

        if mi1 == 0
            continue;  % unmatched stage-1 — retire
        end

        s1     = tracker.currState(buf1_tracks{ti1});
        m2     = unmatched_from_2{mi1};
        x_pos1 = s1.state(s1.state_space.pos_idx);

        if norm(x_pos1) < 1.0 || isempty(msmt_model.least_square_fun)
            continue;
        end

        dt = m2.time - s1.time;
        if dt <= 0, continue; end

        % Second-point LS seeded from stage-1 position
        [x_pos2, ~] = msmt_model.least_square_fun(m2.zeta, x_pos1);
        if any(isnan(x_pos2(:))) || norm(x_pos2) < 1.0 || norm(x_pos2) > 5e6
            continue;
        end

        ss = motion_model.state_space;
        x2 = zeros(ss.num_states, 1);
        x2(ss.pos_idx) = x_pos2;
        % Velocity block: rough finite-difference estimate for gating only
        vel_rough = (x_pos2 - x_pos1) / dt;
        if ss.has_vel
            x2(ss.vel_idx) = vel_rough;
        end
        P2 = build_initial_covariance_two_point(ss, x_pos2, vel_rough, m2.time, dt, ...
                                                msmt_model, target_max_velocity, ...
                                                target_max_acceleration);
        s2   = tracker.makeState(ss, m2.time, x2, P2);
        trk2 = tracker.appendTrack(buf1_tracks{ti1}, s2, false);  % inherits track_id
        new_buf2_msmts{end+1}  = m2;    %#ok<AGROW>
        new_buf2_tracks{end+1} = trk2;  %#ok<AGROW>

        if verbose
            fprintf('  [3PT stage2] trk_id=%d  s1=[%.1f %.1f %.1f] km -> s2=[%.1f %.1f %.1f] km\n', ...
                trk2.track_id, x_pos1(1)/1e3, x_pos1(2)/1e3, x_pos1(3)/1e3, ...
                x_pos2(1)/1e3, x_pos2(2)/1e3, x_pos2(3)/1e3);
        end
    end

    unmatched_from_1 = unmatched_from_2(unmatched_msmt_idx1);
end
% Stage-1 tracks are ALL retired one-shot
new_buf1_msmts  = {};
new_buf1_tracks = {};

%% Step 3: buffer remaining measurements as new stage-1 tracks
for jj = 1:numel(unmatched_from_1)
    [t_buf, ok] = build_single_point_track(unmatched_from_1{jj}, motion_model, msmt_model, ...
                                           target_max_velocity, target_max_acceleration, ...
                                           next_track_id);
    if ~ok
        continue;
    end
    next_track_id = next_track_id + 1;
    new_buf1_msmts{end+1}  = unmatched_from_1{jj};  %#ok<AGROW>
    new_buf1_tracks{end+1} = t_buf;                  %#ok<AGROW>
end

end  % initiate_three_point


%% ==========================================================================
%% Shared local helpers
%% ==========================================================================

function [t_buf, ok] = build_single_point_track(m, motion_model, msmt_model, ...
                                                  target_max_velocity, ...
                                                  target_max_acceleration, track_id)
% Attempt to build a single-point Track from a measurement via LS.
% ok = false if LS failed (position near origin); t_buf is still returned
% but should not be used when ok = false.

ss    = motion_model.state_space;
n_pos = numel(ss.pos_idx);
x_buf = zeros(ss.num_states, 1);
P_buf = 1e6 * eye(ss.num_states);
ok    = false;

if ~isempty(msmt_model.least_square_fun)
    try
        [x_pos_est, ~] = msmt_model.least_square_fun(m.zeta, zeros(n_pos, 1));
        if any(isnan(x_pos_est(:)))
            error('tracker:initiateTracks:lsNaN', 'LS returned NaN');
        end
        % z-flip retry to escape the mirror-image local minimum from
        % near-coplanar sensor arrays
        if n_pos >= 3 && x_pos_est(3) < 0
            x_retry = zeros(n_pos, 1);
            x_retry(3) = -x_pos_est(3);
            [x_retry2, ~] = msmt_model.least_square_fun(m.zeta, x_retry);
            if x_retry2(3) >= 0
                x_pos_est = x_retry2;
            end
        end
        if norm(x_pos_est) >= 1.0 && norm(x_pos_est) < 5e5
            x_buf(ss.pos_idx) = x_pos_est;
            if ~isempty(msmt_model.crlb_fun)
                crlb_pos = msmt_model.crlb_fun(x_pos_est);
                if ~any(isnan(crlb_pos(:)))
                    P_buf(ss.pos_idx, ss.pos_idx) = crlb_pos;
                end
            end
            P_buf = apply_covar_caps(P_buf, ss, target_max_velocity, target_max_acceleration);
            ok = true;
        end
    catch
        % LS failed or diverged — ok stays false
    end
end

s_buf = tracker.makeState(ss, m.time, x_buf, P_buf);
t_buf = tracker.makeTrack(s_buf, track_id);
t_buf.max_velocity     = target_max_velocity;
t_buf.max_acceleration = target_max_acceleration;

end  % build_single_point_track


function P = apply_covar_caps(P, ss, target_max_velocity, target_max_acceleration)
% Scale velocity/acceleration covariance blocks so no diagonal exceeds bound².

if ~isempty(target_max_velocity) && ss.has_vel && ~isempty(ss.vel_idx)
    vel_block = P(ss.vel_idx, ss.vel_idx);
    max_diag  = max(diag(vel_block));
    if max_diag > target_max_velocity^2
        P(ss.vel_idx, ss.vel_idx) = vel_block * (target_max_velocity^2 / max_diag);
    end
end
if ~isempty(target_max_acceleration) && isfield(ss, 'accel_idx') && ~isempty(ss.accel_idx)
    acc_block = P(ss.accel_idx, ss.accel_idx);
    max_diag  = max(diag(acc_block));
    if max_diag > target_max_acceleration^2
        P(ss.accel_idx, ss.accel_idx) = acc_block * (target_max_acceleration^2 / max_diag);
    end
end

end  % apply_covar_caps


function P = build_initial_covariance_two_point(ss, x_pos, vel_est, t, dt, msmt_model, ...
                                                 target_max_velocity, target_max_acceleration)
% Two-point covariance:
%   P_pos   = 10 * CRLB(x_pos)
%   P_vel   = P_pos / dt^2
%   P_accel = P_vel / dt^2

pos_covar_multiplier = 10.0;

if ~isempty(msmt_model.crlb_fun)
    crlb_raw = msmt_model.crlb_fun(x_pos);
    if any(isnan(crlb_raw(:)))
        P = 1e6 * eye(ss.num_states);
        return;
    end
    crlb_pos = pos_covar_multiplier * crlb_raw;
else
    x_temp = zeros(ss.num_states, 1);
    x_temp(ss.pos_idx) = x_pos;
    if ss.has_vel && ~isempty(ss.vel_idx)
        x_temp(ss.vel_idx) = vel_est;
    end
    s_temp = tracker.makeState(ss, t, x_temp, eye(ss.num_states));
    H = msmt_model.h_fun(s_temp);
    if ndims(H) == 3, H = H(:,:,1); end
    H_pos = H(:, ss.pos_idx);
    FIM   = H_pos' * (msmt_model.R \ H_pos);
    if rcond(FIM) < 1e-12
        P = 1e6 * eye(ss.num_states);
        return;
    end
    crlb_pos = pos_covar_multiplier * inv(FIM);  %#ok<MINV>
end

crlb_vel   = crlb_pos / (dt^2);
crlb_accel = crlb_vel  / (dt^2);

P = zeros(ss.num_states);
P(ss.pos_idx, ss.pos_idx) = crlb_pos;
if ss.has_vel && ~isempty(ss.vel_idx)
    P = set_cov_block(P, ss.vel_idx, crlb_vel, target_max_velocity);
end
if isfield(ss, 'accel_idx') && ~isempty(ss.accel_idx)
    P = set_cov_block(P, ss.accel_idx, crlb_accel, target_max_acceleration);
end

end  % build_initial_covariance_two_point


function P = build_initial_covariance_three_point(ss, x_pos1, x_pos2, x_pos3, dt, ...
                                                   msmt_model, target_max_velocity, ...
                                                   target_max_acceleration)
% Three-point covariance (central / second finite-difference propagation):
%   pos_var   = 10 * mean(CRLB at each of the three stages)
%   P_vel     = pos_var / (2 * dt^2)   [central difference]
%   P_accel   = 6 * pos_var / dt^4     [second difference]

pos_covar_multiplier = 10.0;

if isempty(msmt_model.crlb_fun)
    P = 1e6 * eye(ss.num_states);
    return;
end

crlb1 = msmt_model.crlb_fun(x_pos1);
crlb2 = msmt_model.crlb_fun(x_pos2);
crlb3 = msmt_model.crlb_fun(x_pos3);
if any(isnan(crlb1(:))) || any(isnan(crlb2(:))) || any(isnan(crlb3(:)))
    P = 1e6 * eye(ss.num_states);
    return;
end

pos_var    = pos_covar_multiplier * (crlb1 + crlb2 + crlb3) / 3;
crlb_vel   = pos_var / (2.0 * dt^2);
crlb_accel = 6.0 * pos_var / (dt^4);

P = zeros(ss.num_states);
P(ss.pos_idx, ss.pos_idx) = pos_var;
if ss.has_vel && ~isempty(ss.vel_idx)
    P = set_cov_block(P, ss.vel_idx, crlb_vel, target_max_velocity);
end
if isfield(ss, 'accel_idx') && ~isempty(ss.accel_idx)
    P = set_cov_block(P, ss.accel_idx, crlb_accel, target_max_acceleration);
end

end  % build_initial_covariance_three_point


function P = set_cov_block(P, idx, cov_block, max_val)
% Write cov_block into P(idx,idx), capping so max diagonal <= max_val^2.
if ~isempty(max_val)
    max_diag = max(diag(cov_block));
    if max_diag > max_val^2
        cov_block = cov_block * (max_val^2 / max_diag);
    end
end
P(idx, idx) = cov_block;
end  % set_cov_block
