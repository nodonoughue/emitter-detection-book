function [ve, vn, vu] = ecef2enu_vel(vx, vy, vz, lat0, lon0, angle_units)
% [ve, vn, vu] = ecef2enu_vel(vx, vy, vz, lat0, lon0, angle_units, dist_units)
%
% Convert a set of ECEF velocities to ENU (east, north, up),
% relative to a reference Lat/Lon point.  ECEF inputs must be 
% broadcastable to a common size (all non-singleton dimensions must match).
%
% Similar to ecef2enu, except that the origin is not translated; only
% rotation occurs, since the velocity vector is not a pointing vector from
% the ENU origin to the target, rather it is a pointing vector for how
% the target is moving.
%
% The reference lat/lon must be scalar.
%
% The optional input angle_units is used to specify the units for lat/lon 
% (either radians or degrees)
%
% INPUTS:
%   vx              Vector or matrix of ECEF velocity x-components
%   vy              Vector or matrix of ECEF velocity y-components
%   vz              Vector or matrix of ECEF velocity z-components
%   lat0            Reference latitude
%   lon0            Reference longitude
%   alt0            Reference altitude
%   angle_units     Units for lat and lon [Default='deg']
%
% OUTPUTS:
%   ve              East component of velocity (m)
%   vn              North component of velocity(m)
%   vu              Up component of velocity (m)
%
% Nicholas O'Donoughue
% 24 June 2021

if nargin < 7 || isempty(angle_units)
    angle_units = 'deg';
end

% Compute some sin/cos terms
angle_mult = unitsratio('rad', angle_units);
lat_rad = angle_mult * lat0;
lon_rad = angle_mult * lon0;

cos_lat = cos(lat_rad);
cos_lon = cos(lon_rad);
sin_lat = sin(lat_rad);
sin_lon = sin(lon_rad);

% Compute Compensation term (t)
vt = (cos_lon .* vx) + (sin_lon .* vy);

% Compute e, n, and u
ve = -(sin_lon .* vx) + (cos_lon .* vy);
vn = -(sin_lat .* vt)  + (cos_lat .* vz);
vu =  (cos_lat .* vt)  + (sin_lat .* vz);