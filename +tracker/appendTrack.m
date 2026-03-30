function trk = appendTrack(trk, state, missed_detection)
% appendTrack  Append a new state to an existing track.
%
% trk = appendTrack(trk, state)
% trk = appendTrack(trk, state, missed_detection)
%
% INPUTS
%   trk               Track struct from makeTrack
%   state             State struct from makeState to append
%   missed_detection  Logical flag.  If true, increments
%                     num_missed_detections and does NOT increment
%                     num_updates.  If false (default), resets
%                     num_missed_detections to 0 and increments
%                     num_updates.
%
% OUTPUTS
%   trk   Updated track struct
%
% Nicholas O'Donoughue
% June 2025

if nargin < 3 || isempty(missed_detection)
    missed_detection = false;
end

trk.states{end+1} = state;

if missed_detection
    trk.num_missed_detections = trk.num_missed_detections + 1;
else
    trk.num_missed_detections = 0;
    trk.num_updates = trk.num_updates + 1;
end
