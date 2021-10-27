function [vx, vy, vz] = enu2ecef_vel(ve, vn, vu, lat0, lon0, alt0, angle_units)
% [vx, vy, vz] = enu2ecef_vel(ve, vn, vu, lat0, lon0, alt0, angle_units)
%
% Convert local ENU velocity vector to ECEF, given the reference LLA position
% for the local ENU coordinates.
%
% Similar to enu2ecef, except that this function does not translate to
% account for migration of the origin from lat0, lon0 to the center of the
% Earth, since the input vectors are taken as referenced from target
% position.  Instead, only the rotation operation is conducted.
%
% INPUTS:
%   ve              East component of velocity(m)
%   vn              North component of velocity (m)
%   vu              Up component of velocity (m)
%   lat0            Reference latitude
%   lon0            Reference longitude
%   alt0            Reference altitude
%   angle_units     Units for lat and lon [Default='deg']
%
% OUTPUTS:
%   vx              ECEF velocity x-component (m)
%   vy              ECEF velocity y-component(m)
%   vz              ECEF velocity z-component(m)
%
% Nicholas O'Donoughue
% 9 June 2021

if nargin < 7 || isempty(angle_units)
    angle_units = 'deg';
end

% Precompute Trig Functions of Reference Lat/Lon
lat_rad = lat0 * unitsratio('rad', angle_units);
lon_rad = lon0 * unitsratio('rad', angle_units);

sin_lat = sin(lat_rad);
sin_lon = sin(lon_rad);
cos_lat = cos(lat_rad);
cos_lon = cos(lon_rad);

% Convert from ENU to dx/dy/dz, using reference coordinates
vt = -sin_lat.*vn + cos_lat.*vu;

vx = ((-sin_lon .* ve) + (cos_lon .* vt));
vy = ((cos_lon .* ve) + (sin_lon .* vt));
vz = ((cos_lat .* vn) + (sin_lat .* vu));
