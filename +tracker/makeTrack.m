function trk = makeTrack(initial_state, track_id, motion_model)
% makeTrack  Create a new tracker track struct.
%
% trk = makeTrack(initial_state)
% trk = makeTrack(initial_state, track_id)
% trk = makeTrack(initial_state, track_id, motion_model)
%
% INPUTS
%   initial_state  State struct from makeState (the seed state for the track)
%   track_id       Optional identifier (string or integer).  Default: ''
%   motion_model   Optional motion model struct from makeMotionModel.
%                  When set, makeHypothesis and associateTracks infer the
%                  model from the track rather than requiring it as an
%                  explicit argument.  Default: []
%
% OUTPUTS
%   trk   Struct with fields:
%           states                – 1-D cell array of State structs
%           track_id              – identifier (string or integer)
%           num_missed_detections – consecutive missed-detection counter
%           num_updates           – total successful updates received
%           max_velocity          – optional speed cap [m/s]; [] = unconstrained
%           max_acceleration      – optional acceleration cap [m/s²]; [] = unconstrained
%           motion_model          – motion model struct (or [] if not set)
%
% Nicholas O'Donoughue
% June 2025

if nargin < 2 || isempty(track_id)
    track_id = '';
end
if nargin < 3
    motion_model = [];
end

trk = struct('states',                {{initial_state}}, ...
             'track_id',              track_id, ...
             'num_missed_detections', 0, ...
             'num_updates',           1, ...
             'max_velocity',          [], ...
             'max_acceleration',      [], ...
             'motion_model',          motion_model);
