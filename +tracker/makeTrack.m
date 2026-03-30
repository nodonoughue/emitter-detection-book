function trk = makeTrack(initial_state, track_id)
% makeTrack  Create a new tracker track struct.
%
% trk = makeTrack(initial_state)
% trk = makeTrack(initial_state, track_id)
%
% INPUTS
%   initial_state  State struct from makeState (the seed state for the track)
%   track_id       Optional identifier (string or integer).  Default: ''
%
% OUTPUTS
%   trk   Struct with fields:
%           states                – 1-D cell array of State structs
%           track_id              – identifier (string or integer)
%           num_missed_detections – consecutive missed-detection counter
%           num_updates           – total successful updates received
%
% Nicholas O'Donoughue
% June 2025

if nargin < 2 || isempty(track_id)
    track_id = '';
end

trk = struct('states',                {{initial_state}}, ...
             'track_id',              track_id, ...
             'num_missed_detections', 0, ...
             'num_updates',           1);
