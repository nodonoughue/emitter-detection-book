function [x, y, z] = aer2ecef(a, e, r, lat0, lon0, alt0, angle_units, dist_units)
% [x, y, z] = aer2ecef(a, e, r, angle_units, dist_units)
%
% Convert spherical AER coordinates to ECEF, where azimuth is
% in degrees East from North (as opposed to the typical degrees +y from +x)
% and elevation is degrees above the horizontal.
%
% INPUTS:
%   a               Azimuth (deg E from N)
%   e               Elevation (deg above horizontal)
%   r               Range (meters)
%   lat0            Reference latitude
%   lon0            Reference longitude
%   alt0            Reference altitude
%   angle_units     Units for lat and lon [Default='deg']
%   dist_units      Units for alt [Default='m']
%
% OUTPUTS:
%   x               ECEF x-coordinate (m)
%   y               ECEF y-coordinate (m)
%   z               ECEF z-coordinate (m)
%
% Nicholas O'Donoughue
% 9 June 2021

if nargin < 7 || isempty(angle_units)
    angle_units = 'deg';
end

if nargin < 8 || isempty(dist_units)
    dist_units = 'm';
end

% Convert to ENU
[e, n, u] = utils.aer2enu(a,e,r,angle_units,dist_units);

% Ensure alt0 is meters, to match ENU distance units
alt0 = alt0 * unitsratio('m', dist_units);

% Convert from ENU to ECEF
[x, y, z] = utils.enu2ecef(e, n, u, lat0, lon0,alt0,angle_units, 'm');