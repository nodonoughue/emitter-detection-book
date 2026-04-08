function s = constrainMotion(s, max_velocity, max_acceleration)
% constrainMotion  Clip velocity and/or acceleration so that neither exceeds
%                  the given magnitude bound.  Direction is preserved; only
%                  the magnitude is scaled down.
%
% s = constrainMotion(s, max_velocity, max_acceleration)
%
% INPUTS
%   s                State struct from makeState
%   max_velocity     Maximum speed [m/s], or [] to skip the velocity check
%   max_acceleration Maximum acceleration magnitude [m/s²], or [] to skip
%
% OUTPUTS
%   s   State struct with velocity and/or acceleration clipped in-place
%
% Nicholas O'Donoughue
% June 2025

ss = s.state_space;

if ~isempty(max_velocity) && ss.has_vel && ~isempty(ss.vel_idx)
    v     = s.state(ss.vel_idx);
    speed = norm(v);
    if speed > max_velocity
        s.state(ss.vel_idx) = v * (max_velocity / speed);
    end
end

if ~isempty(max_acceleration) && isfield(ss, 'accel_idx') && ~isempty(ss.accel_idx)
    a   = s.state(ss.accel_idx);
    mag = norm(a);
    if mag > max_acceleration
        s.state(ss.accel_idx) = a * (max_acceleration / mag);
    end
end
