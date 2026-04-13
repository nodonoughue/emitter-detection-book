function tracker_state = makeTrackerState(motion_model, msmt_model, varargin)
% makeTrackerState  Create the mutable state struct used by runTrackerStep.
%
% tracker_state = makeTrackerState(motion_model, msmt_model)
% tracker_state = makeTrackerState(motion_model, msmt_model, Name, Value, ...)
%
% INPUTS
%   motion_model  Motion model struct from makeMotionModel
%   msmt_model    Measurement model struct from makeMeasurementModel
%
%   Name-Value pairs (all optional):
%     'gate_probability'     Chi-square gate probability (default: 0.99)
%     'max_missed'           Max consecutive missed detections (default: 3)
%     'num_hits'             M-of-N promoter: required hit count (default: 3)
%     'num_chances'          M-of-N promoter: total chances before discard (default: 5)
%     'assoc_type'           'nn' or 'gnn' (default: 'gnn')
%     'keep_all_tracks'      Logical, keep deleted tracks for later analysis (default: true)
%     'target_max_velocity'     Max target speed [m/s] for covariance caps and state
%                               clipping after EKF updates.  [] = unconstrained (default: [])
%     'target_max_acceleration' Max target acceleration [m/s²] for covariance caps and
%                               state clipping.  Only effective for CA/CJ state spaces.
%                               [] = unconstrained (default: [])
%
% OUTPUTS
%   tracker_state  Struct with fields:
%                    config           – struct of configuration parameters
%                    firm_tracks      – {} (initially empty)
%                    tentative_tracks – {} (initially empty)
%                    deleted_tracks   – {} (initially empty)
%                    buffer_msmts     – {} (initiator measurement buffer)
%                    buffer_tracks    – {} (initiator track buffer)
%                    next_track_id    – 1
%
% Nicholas O'Donoughue
% June 2025

% --- Parse optional parameters ---
p = inputParser;
addParameter(p, 'gate_probability',      0.99);
addParameter(p, 'max_missed',            3);
addParameter(p, 'num_hits',              3);
addParameter(p, 'num_chances',           5);
addParameter(p, 'assoc_type',            'gnn');
addParameter(p, 'detection_probability', 1.0);
addParameter(p, 'keep_all_tracks',       true);
addParameter(p, 'target_max_velocity',     []);
addParameter(p, 'target_max_acceleration', []);
addParameter(p, 'verbose',                 false);
parse(p, varargin{:});
opt = p.Results;

config = struct('motion_model',            motion_model, ...
                'msmt_model',              msmt_model, ...
                'gate_probability',        opt.gate_probability, ...
                'max_missed',              opt.max_missed, ...
                'num_hits',                opt.num_hits, ...
                'num_chances',             opt.num_chances, ...
                'assoc_type',              opt.assoc_type, ...
                'detection_probability',   opt.detection_probability, ...
                'keep_all_tracks',         opt.keep_all_tracks, ...
                'target_max_velocity',     opt.target_max_velocity, ...
                'target_max_acceleration', opt.target_max_acceleration, ...
                'verbose',                 opt.verbose);

tracker_state = struct('config',           config, ...
                       'firm_tracks',      {{}}, ...
                       'tentative_tracks', {{}}, ...
                       'deleted_tracks',   {{}}, ...
                       'failed_tracks',    {{}}, ...
                       'buffer_msmts',     {{}}, ...
                       'buffer_tracks',    {{}}, ...
                       'next_track_id',    1);
