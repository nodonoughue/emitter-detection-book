function [x, y, z] = enu2ecef(e, n, u, lat0, lon0, alt0, angle_units, dist_units)
% [x, y, z] = enu2ecef(e, n, u, lat0, lon0, alt0, angle_units, dist_units)
%
% Convert local ENU coordinates to ECEF, given the reference LLA position
% for the local ENU coordinates.
%
% INPUTS:
%   e               East (m)
%   n               North (m)
%   u               Up (m)
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

% Precompute Trig Functions of Reference Lat/Lon
lat_rad = lat0 * unitsratio('rad', angle_units);
lon_rad = lon0 * unitsratio('rad', angle_units);

sin_lat = sin(lat_rad);
sin_lon = sin(lon_rad);
cos_lat = cos(lat_rad);
cos_lon = cos(lon_rad);

% Convert from ENU to dx/dy/dz, using reference coordinates
t = -sin_lat.*n + cos_lat.*u;

dx = ((-sin_lon .* e) + (cos_lon .* t)) * unitsratio('m', dist_units);
dy = ((cos_lon .* e) + (sin_lon .* t)) * unitsratio('m', dist_units);
dz = ((cos_lat .* n) + (sin_lat .* u)) * unitsratio('m', dist_units);

% Convert Ref LLA to ECEF; output in meters
[x0, y0, z0] = utils.lla2ecef(lat0, lon0, alt0, angle_units, dist_units);

% Combine
x = x0 + dx;
y = y0 + dy;
z = z0 + dz;