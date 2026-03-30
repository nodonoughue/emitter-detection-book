function s = currState(trk)
% currState  Return the most recent state from a track struct.
%
% s = currState(trk)
%
% INPUTS
%   trk   Track struct from makeTrack / appendTrack
%
% OUTPUTS
%   s     The last State struct in trk.states
%
% Nicholas O'Donoughue
% June 2025

s = trk.states{end};
