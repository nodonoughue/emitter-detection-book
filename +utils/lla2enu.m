function [e, n, u] = lla2enu(lat, lon, alt, lat0, lon0, alt0, angle_units, dist_units)
% [x, y, z] = lla2ecef(lat, lon, alt, lat0, lon0, alt0, angle_units, dist_units)
%
% Convert a set of Lat, Lon, Alt coordinates to ENU (east, north, up),
% relative to a reference Lat/Lon point.  Lat, Lon, and Alt inputs must be 
% broadcastable to a common size (all non-singleton dimensions must match).
%
% The optional inputs angle_units and dist_units are used to specify the
% units for lat/lon (either radians or degrees), and alt (any valid length
% unit).  See 'help unitsratio' for information on valid specifications.
%
% INPUTS:
%   lat             Vector or matrix of latitudes
%   lon             Vector or matrix of longitudes
%   alt             Vector or matrix of altitudes
%   lat0            Reference latitude
%   lon0            Reference longitude
%   alt0            Reference altitude
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

if nargin < 7 || isempty(angle_units)
    angle_units = 'deg';
end

if nargin < 8 || isempty(dist_units)
    dist_units = 'm';
end

% First, we convert the points to ECEF. The result will be in meters
[x, y, z] = utils.lla2ecef(lat, lon, alt, angle_units, dist_units);

% Make sure ref altitude is in meters (to match ECEF distance units)
alt0 = alt0 * unitsratio('m', dist_units);

% Then, we leverage ecef2enu to compute
[e, n, u] = utils.ecef2enu(x, y, z, lat0, lon0, alt0, angle_units, 'm');
