function [e, n, u] = ecef2enu(x, y, z, lat0, lon0, alt0, angle_units, dist_units)
% [e, n, u] = ecef2enu(x, y, z, lat0, lon0, alt0, angle_units, dist_units)
%
% Convert a set of ECEF coordinates to ENU (east, north, up),
% relative to a reference Lat/Lon point.  ECEF inputs must be 
% broadcastable to a common size (all non-singleton dimensions must match).
%
% The reference lat/lon must be scalar.
%
% The optional inputs angle_units and dist_units are used to specify the
% units for lat/lon (either radians or degrees), and alt (any valid length
% unit).  See 'help unitsratio' for information on valid specifications.
%
% INPUTS:
%   x               Vector or matrix of ECEF x-coordinates
%   y               Vector or matrix of ECEF y-coordinates
%   z               Vector or matrix of ECEF z-coordinates
%   lat0            Reference latitude
%   lon0            Reference longitude
%   alt0            Reference altitude
%   angle_units     Units for lat and lon [Default='deg']
%   dist_units      Units for ECEF coordinates [Default='m']
%
% OUTPUTS:
%   e               East (m)
%   n               North (m)
%   u               Up (m)
%
% Nicholas O'Donoughue
% 9 June 2021

if nargin < 7 || isempty(angle_units)
    angle_units = 'deg';
end

if nargin < 8 || isempty(dist_units)
    dist_units = 'm';
end

% First, we convert the reference point to ECEF. The result will be in
% meters.
[x0, y0, z0] = utils.lla2ecef(lat0, lon0, alt0, angle_units, dist_units);

% Take the difference
dx = x - x0;
dy = y - y0;
dz = z - z0;

% Compute some sin/cos terms
angle_mult = unitsratio('rad', angle_units);
lat_rad = angle_mult * lat0;
lon_rad = angle_mult * lon0;

cos_lat = cos(lat_rad);
cos_lon = cos(lon_rad);
sin_lat = sin(lat_rad);
sin_lon = sin(lon_rad);

% Compute Compensation term (t)
t = (cos_lon .* dx) + (sin_lon .* dy);

% Compute e, n, and u
e = -(sin_lon .* dx) + (cos_lon .* dy);
n = -(sin_lat .* t)  + (cos_lat .* dz);
u =  (cos_lat .* t)  + (sin_lat .* dz);