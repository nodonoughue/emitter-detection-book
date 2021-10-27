function [lat, lon, alt] = enu2lla(e, n, u, lat0, lon0, alt0, angle_units, dist_units)
% [lat, lon, alt] = enu2lla(e, n, u, lat0, lon0, alt0, angle_units, dist_units)
%
% Convert local ENU coordinates to LLA, given the reference LLA position
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
%   lat             Latitude (degrees North)
%   lon             Longitude (degrees East)
%   alt             Altitude (meters)
%
% Nicholas O'Donoughue
% 9 June 2021

if nargin < 7 || isempty(angle_units)
    angle_units = 'deg';
end

if nargin < 8 || isempty(dist_units)
    dist_units = 'm';
end

% Convert to global ECEF coordinates
[x, y, z] = utils.enu2ecef(e, n, u, lat0, lon0, alt0, angle_units, dist_units);

% Convert from ECEF to LLA
[lat, lon, alt] = utils.ecef2lla(x, y, z, 'm');