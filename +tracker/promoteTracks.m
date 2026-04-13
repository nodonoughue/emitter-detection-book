function [to_promote, to_delete, to_keep] = promoteTracks(tracks, num_hits, num_chances, verbose)
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

if nargin < 4, verbose = false; end

to_promote = {};
to_delete  = {};
to_keep    = {};

for ii = 1:numel(tracks)
    t = tracks{ii};
    if t.num_updates >= num_hits
        to_promote{end+1} = t;  %#ok<AGROW>
        if verbose
            s   = tracker.currState(t);
            ss  = s.state_space;
            pos = s.state(ss.pos_idx);
            vel = [];
            if ss.has_vel && ~isempty(ss.vel_idx)
                vel = s.state(ss.vel_idx);
            end
            fprintf('  [PROMOTE] track_id=%d  updates=%d/%d  pos=[%.1f, %.1f, %.1f] km', ...
                t.track_id, t.num_updates, numel(t.states), ...
                pos(1)/1e3, pos(2)/1e3, pos(3)/1e3);
            if ~isempty(vel)
                fprintf('  vel=[%.1f, %.1f, %.1f] m/s (speed=%.1f)', ...
                    vel(1), vel(2), vel(3), norm(vel));
            end
            fprintf('\n');
        end
    elseif numel(t.states) >= num_chances
        to_delete{end+1} = t;   %#ok<AGROW>
        if verbose
            s   = tracker.currState(t);
            pos = s.state(s.state_space.pos_idx);
            fprintf('  [DROP]    track_id=%d  updates=%d/%d  pos=[%.1f, %.1f, %.1f] km\n', ...
                t.track_id, t.num_updates, numel(t.states), ...
                pos(1)/1e3, pos(2)/1e3, pos(3)/1e3);
        end
    else
        to_keep{end+1} = t;     %#ok<AGROW>
    end
end
