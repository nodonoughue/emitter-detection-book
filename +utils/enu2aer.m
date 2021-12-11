function [a, e, r] = enu2aer(e, n, u, dist_units)
% [a, e, r] = enu2aer(e, n, u, dist_units)
%
% Convert cartesian ENU coordinates to spherical AER, where azimuth is
% in degrees East from North (as opposed to the typical degrees +y from +x)
% and elevation is degrees above the horizontal.
%
% INPUTS:
%   e               East (m)
%   n               North (m)
%   u               Up (m)
%   dist_units      Units for alt [Default='m']
%
% OUTPUTS:
%   a               Azimuth (deg E from N)
%   e               Elevation (deg above horizontal)
%   r               Range (meters)
%
% Nicholas O'Donoughue
% 9 June 2021

if nargin < 4 || isempty(dist_units)
    dist_units = 'm';
end

% Compute Ground and Slant Ranges
ground_range = sqrt(e.^2 + n.^2);
slant_range = sqrt(ground_range.^2 + u.^2);

% Az/El Angles, in Radians
el_angle_rad = atan2(u, ground_range);
az_angle_rad = atan2(e, n);

% Outputs
a = az_angle_rad * 180/pi;
e = el_angle_rad * 180/pi;
r = slant_range * unitsratio('m', dist_units);