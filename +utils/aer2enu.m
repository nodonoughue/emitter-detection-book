function [e, n, u] = aer2enu(a, e, r, angle_units, dist_units)
% [e, n, u] = aer2enu(a, e, r, angle_units, dist_units)
%
% Convert spherical AER coordinates to cartesian ENU, where azimuth is
% in degrees East from North (as opposed to the typical degrees +y from +x)
% and elevation is degrees above the horizontal.
%
% INPUTS:
%   a               Azimuth (deg E from N)
%   e               Elevation (deg above horizontal)
%   r               Range (meters)
%   angle_units     Units for lat and lon [Default='deg']
%   dist_units      Units for alt [Default='m']
%
% OUTPUTS:
%   e               East (m)
%   n               North (m)
%   u               Up (m)
%
% Nicholas O'Donoughue
% 9 June 2021

if nargin < 4 || isempty(angle_units)
    angle_units = 'deg';
end

if nargin < 5 || isempty(dist_units)
    dist_units = 'm';
end

% Fix input units
az_angle_rad = a .* unitsratio('rad', angle_units);
el_angle_rad = e .* unitsratio('rad', angle_units);
range_m = r .* unitsratio('m', dist_units);

% Compute Ground Range
ground_range = range_m .* cos(el_angle_rad);

% Outputs
e = ground_range .* sin(az_angle_rad);
n = ground_range .* cos(az_angle_rad);
u = range_m .* sin(el_angle_rad);