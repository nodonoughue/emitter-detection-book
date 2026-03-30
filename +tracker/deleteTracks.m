function [active, deleted] = deleteTracks(tracks, max_missed)
% deleteTracks  Missed-detection deleter: remove tracks with too many consecutive
%               missed detections.
%
% [active, deleted] = deleteTracks(tracks, max_missed)
%
% INPUTS
%   tracks      Cell array of Track structs (firm or tentative)
%   max_missed  Maximum number of consecutive missed detections allowed;
%               tracks with num_missed_detections > max_missed are deleted
%
% OUTPUTS
%   active   Cell array of tracks that survive (num_missed_detections <= max_missed)
%   deleted  Cell array of tracks that were deleted
%
% Nicholas O'Donoughue
% June 2025

active  = {};
deleted = {};

for ii = 1:numel(tracks)
    if tracks{ii}.num_missed_detections > max_missed
        deleted{end+1} = tracks{ii};  %#ok<AGROW>
    else
        active{end+1}  = tracks{ii};  %#ok<AGROW>
    end
end
