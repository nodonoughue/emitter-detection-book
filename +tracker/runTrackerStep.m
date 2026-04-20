function tracker_state = runTrackerStep(tracker_state, measurements, curr_time)
% runTrackerStep  Process one scan of measurements through the full tracker
%                 pipeline: associate → promote → initiate → delete.
%
% tracker_state = runTrackerStep(tracker_state, measurements, curr_time)
%
% INPUTS
%   tracker_state  Struct produced by makeTrackerState (updated in place)
%   measurements   Cell array of Measurement structs for the current scan.
%                  Each struct must have a non-empty .msmt_model field
%                  (i.e. created via makeMeasurement(msmt_model, state, time)
%                  or makeMeasurement(time, zeta, msmt_model)).
%   curr_time      Timestamp of the current scan [s]; inferred from
%                  measurements{1}.time if empty
%
% OUTPUTS
%   tracker_state  Updated tracker state struct
%
% See also makeTrackerState, tracker.associateTracks, tracker.promoteTracks,
%          tracker.initiateTracks, tracker.deleteTracks
%
% Nicholas O'Donoughue
% June 2025

if nargin < 3 || isempty(curr_time)
    if ~isempty(measurements)
        curr_time = measurements{1}.time;
    else
        curr_time = [];
    end
end

% Validate that every measurement carries a measurement model
for ii = 1:numel(measurements)
    if isempty(measurements{ii}.msmt_model)
        error('tracker:runTrackerStep:missingMsmtModel', ...
              'measurements{%d}.msmt_model is empty. All measurements must be created with a msmt_model.', ii);
    end
end

cfg = tracker_state.config;

%% --- 1. Associate measurements with firm tracks --------------------------
[tracker_state, unassoc_msmts] = update_firm_tracks(tracker_state, measurements, curr_time, cfg);

%% --- 2. Associate remaining measurements with tentative tracks -----------
[tracker_state, unassoc_msmts_2] = update_tentative_tracks(tracker_state, unassoc_msmts, curr_time, cfg);

%% --- 3. Initiate new tentative tracks from leftover measurements ---------
% unassoc_msmts_2 is already the set not consumed by firm or tentative tracks
leftover_for_init = unassoc_msmts_2;

verbose = isfield(cfg, 'verbose') && cfg.verbose;
[new_tracks, tracker_state.next_track_id, ...
 tracker_state.buffer_msmts, tracker_state.buffer_tracks, ...
 tracker_state.stage2_msmts, tracker_state.stage2_tracks] = ...
    tracker.initiateTracks(leftover_for_init, curr_time, ...
                           cfg.motion_model, cfg.msmt_model, ...
                           cfg.gate_probability, tracker_state.next_track_id, ...
                           tracker_state.buffer_msmts, tracker_state.buffer_tracks, ...
                           cfg.target_max_velocity, cfg.target_max_acceleration, verbose, ...
                           cfg.init_type, ...
                           tracker_state.stage2_msmts, tracker_state.stage2_tracks);

tracker_state.tentative_tracks = [tracker_state.tentative_tracks, new_tracks];

%% --- 4. Delete firm and tentative tracks with too many missed detections --
[tracker_state.firm_tracks, del_firm] = ...
    tracker.deleteTracks(tracker_state.firm_tracks, cfg.max_missed);

[tracker_state.tentative_tracks, del_tent] = ...
    tracker.deleteTracks(tracker_state.tentative_tracks, cfg.max_missed);

if cfg.keep_all_tracks
    tracker_state.deleted_tracks = [tracker_state.deleted_tracks, del_firm, del_tent];
end

end


%% --------------------------------------------------------------------------
function [tracker_state, unassoc] = update_firm_tracks(tracker_state, measurements, curr_time, cfg)
% Associate + update firm tracks; return unassociated measurements.

unassoc = measurements;
if isempty(tracker_state.firm_tracks) || isempty(measurements)
    if ~isempty(tracker_state.firm_tracks) && ~isempty(curr_time)
        % Coast all firm tracks
        for ii = 1:numel(tracker_state.firm_tracks)
            trk = tracker_state.firm_tracks{ii};
            s = tracker.currState(trk);
            s_pred = tracker.predictState(s, curr_time, cfg.motion_model);
            s_pred = tracker.constrainMotion(s_pred, trk.max_velocity, trk.max_acceleration);
            tracker_state.firm_tracks{ii} = tracker.appendTrack(trk, s_pred, true);
        end
    end
    return;
end

[trk_idx, msmt_idx, unassigned_msmt_idx, pda_states] = ...
    tracker.associateTracks(tracker_state.firm_tracks, measurements, curr_time, ...
                            cfg.motion_model, cfg.msmt_model, cfg.gate_probability, ...
                            cfg.assoc_type, cfg.detection_probability);

% Update assigned tracks
for kk = 1:numel(trk_idx)
    ti = trk_idx(kk);
    mi = msmt_idx(kk);
    s  = tracker.currState(tracker_state.firm_tracks{ti});
    s_pred = tracker.predictState(s, curr_time, cfg.motion_model);

    trk = tracker_state.firm_tracks{ti};
    if ~isempty(pda_states)
        % PDA: fused state already computed by associateTracks
        s_upd = pda_states{ti};
        s_upd = tracker.constrainMotion(s_upd, trk.max_velocity, trk.max_acceleration);
        tracker_state.firm_tracks{ti} = tracker.appendTrack(trk, s_upd, false);
    elseif mi > 0
        s_upd = tracker.ekfUpdate(s_pred, measurements{mi});
        s_upd = tracker.constrainMotion(s_upd, trk.max_velocity, trk.max_acceleration);
        tracker_state.firm_tracks{ti} = tracker.appendTrack(trk, s_upd, false);
    else
        % Missed detection: coast — constrain to prevent unbounded velocity drift
        s_pred = tracker.constrainMotion(s_pred, trk.max_velocity, trk.max_acceleration);
        tracker_state.firm_tracks{ti} = tracker.appendTrack(trk, s_pred, true);
    end
end

unassoc = measurements(unassigned_msmt_idx);
end


%% --------------------------------------------------------------------------
function [tracker_state, unassoc] = update_tentative_tracks(tracker_state, measurements, curr_time, cfg)
% Associate + update tentative tracks; promote / drop as needed.

unassoc = measurements;
if isempty(tracker_state.tentative_tracks)
    return;
end

if isempty(measurements)
    % Coast all tentative tracks (missed detection)
    for ii = 1:numel(tracker_state.tentative_tracks)
        trk = tracker_state.tentative_tracks{ii};
        s = tracker.currState(trk);
        s_pred = tracker.predictState(s, curr_time, cfg.motion_model);
        s_pred = tracker.constrainMotion(s_pred, trk.max_velocity, trk.max_acceleration);
        tracker_state.tentative_tracks{ii} = tracker.appendTrack(trk, s_pred, true);
    end
    % Still run promoter/deleter so expired tentative tracks are culled
    verbose = isfield(cfg, 'verbose') && cfg.verbose;
    [to_promote, to_drop, to_keep] = ...
        tracker.promoteTracks(tracker_state.tentative_tracks, cfg.num_hits, cfg.num_chances, verbose);
    tracker_state.firm_tracks      = [tracker_state.firm_tracks, to_promote];
    tracker_state.tentative_tracks = to_keep;
    if cfg.keep_all_tracks
        tracker_state.failed_tracks = [tracker_state.failed_tracks, to_drop];
    end
    return;
end

[trk_idx, msmt_idx, unassigned_msmt_idx, pda_states] = ...
    tracker.associateTracks(tracker_state.tentative_tracks, measurements, curr_time, ...
                            cfg.motion_model, cfg.msmt_model, cfg.gate_probability, ...
                            cfg.assoc_type, cfg.detection_probability);

% Update tentative tracks
for kk = 1:numel(trk_idx)
    ti = trk_idx(kk);
    mi = msmt_idx(kk);
    s  = tracker.currState(tracker_state.tentative_tracks{ti});
    s_pred = tracker.predictState(s, curr_time, cfg.motion_model);

    trk = tracker_state.tentative_tracks{ti};
    if ~isempty(pda_states)
        % PDA: fused state already computed by associateTracks
        s_upd = pda_states{ti};
        s_upd = tracker.constrainMotion(s_upd, trk.max_velocity, trk.max_acceleration);
        tracker_state.tentative_tracks{ti} = tracker.appendTrack(trk, s_upd, false);
    elseif mi > 0
        s_upd = tracker.ekfUpdate(s_pred, measurements{mi});
        s_upd = tracker.constrainMotion(s_upd, trk.max_velocity, trk.max_acceleration);
        tracker_state.tentative_tracks{ti} = tracker.appendTrack(trk, s_upd, false);
    else
        % Missed detection: coast — constrain to prevent unbounded velocity drift
        s_pred = tracker.constrainMotion(s_pred, trk.max_velocity, trk.max_acceleration);
        tracker_state.tentative_tracks{ti} = tracker.appendTrack(trk, s_pred, true);
    end
end

% Promote / drop tentative tracks
verbose = isfield(cfg, 'verbose') && cfg.verbose;
[to_promote, to_drop, to_keep] = ...
    tracker.promoteTracks(tracker_state.tentative_tracks, cfg.num_hits, cfg.num_chances, verbose);

tracker_state.firm_tracks      = [tracker_state.firm_tracks, to_promote];
tracker_state.tentative_tracks = to_keep;

% Tentative tracks that exhausted their chances without enough hits are
% "failed" tracks — they were never confirmed.  Store separately so that
% deleted_tracks mirrors Python's semantics: only ever-confirmed tracks.
tracker_state.failed_tracks = [tracker_state.failed_tracks, to_drop];

unassoc = measurements(unassigned_msmt_idx);
end
