function [to_promote, to_delete, to_keep] = promoteTracks(tracks, num_hits, num_chances)
% promoteTracks  M-of-N promoter: classify tentative tracks as promote / delete / keep.
%
% [to_promote, to_delete, to_keep] = promoteTracks(tracks, num_hits, num_chances)
%
% INPUTS
%   tracks       Cell array of tentative Track structs
%   num_hits     Minimum num_updates required to promote a track (M)
%   num_chances  Maximum total states allowed before a track is discarded (N)
%
% OUTPUTS
%   to_promote   Cell array of tracks that have >= num_hits updates (promote to firm)
%   to_delete    Cell array of tracks that have >= num_chances states but < num_hits
%                updates (discard)
%   to_keep      Cell array of tracks that have neither threshold yet (keep tentative)
%
% Nicholas O'Donoughue
% June 2025

to_promote = {};
to_delete  = {};
to_keep    = {};

for ii = 1:numel(tracks)
    t = tracks{ii};
    if t.num_updates >= num_hits
        to_promote{end+1} = t;  %#ok<AGROW>
    elseif numel(t.states) >= num_chances
        to_delete{end+1} = t;   %#ok<AGROW>
    else
        to_keep{end+1} = t;     %#ok<AGROW>
    end
end
