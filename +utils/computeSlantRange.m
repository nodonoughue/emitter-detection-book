function R = computeSlantRange(alt1,alt2,el_angle_deg,use_effective_earth)
% R = computeSlantRange(alt1,alt2,el_angle_deg,use_effective_earth)
%
% Computes the slant range between two points given the altitude above the 
% Earth, and the elevation angle relative to ground plane horizontal.
%
% R = sqrt((r1*cos(th))^2 + r2^2 - r1^2) - r1*cos(th)
%
% where
%   r1 = alt1 + Re
%   r2 = alt2 + Re
%   th = elevation angle (above horizon) in degrees at point 1
%
% If the fourth argument is specified and set to true,
% then Re is the 4/3 Earth Radius used for RF propagation paths
% otherwise the true earth radius is used.
%
% INPUTS:
%   alt1                Altitude at the start of the path, in meters
%   alt2                Altitude at the end of the path, in meters
%   el_angle_deg        Elevation angle, at the start of the path, in
%                       degrees above the local horizontal plane
%   use_effective_earth Binary flag [Default=False], specifying whether to
%                       use the 4/3 Earth Radius approximation common to
%                       RF propagation
%
% OUTPUTS:
%   R                   Straight line slant range between the start and 
%                       end point specified, in meters
%
% Nicholas O'Donoughue
% 1 July 2019

if nargin < 4 || isempty(use_effective_earth)
    use_effective_earth = false;
end

if use_effective_earth
    Re = utils.constants.Re;
else
    Re = utils.constants.Re_true;
end

% Compute the two radii
r1 = Re + alt1;
r2 = Re + alt2;
r1c = r1 .* sind(el_angle_deg);

% Compute the slant range
R = sqrt(r1c.^2 + r2.^2 - r1.^2)-r1c;
