function [a, e, r] = lla2aer(lat, lon, alt, lat0, lon0, alt0, angle_units, dist_units)
% [a, e, r] = lla2aer(lat, lon, alt, lat0, lon0, alt0, angle_units, dist_units)
%
% Convert a set of Lat, Lon, Alt coordinates to AER (az, el, range), as
% seen by a reference sensor at lat0, lon0, alt0.  Any non-scalaer LLA
% inputs must be broadcastable to a common shape.
%
% The optional inputs angle_units and dist_units are used to specify the
% units for lat/lon (either radians or degrees), and alt (any valid length
% unit).  See 'help unitsratio' for information on valid specifications.
%
% INPUTS:
%   lat             Vector or matrix of latitudes
%   lon             Vector or matrix of longitudes
%   alt             Vector or matrix of altitudes
%   lat0            Vector or matrix of reference latitudes
%   lon0            Vector or matrix of reference longitudes
%   alt0            Vector or matrix of reference altitudes
%   angle_units     Units for lat and lon [Default='deg']
%   dist_units      Units for alt [Default='m']
%
% OUTPUTS:
%   a               Azimuth (deg E from N)
%   e               Elevation (deg above horizontal)
%   r               Range (meters)
%
% Nicholas O'Donoughue
% 9 June 2021

if nargin < 7 || isempty(angle_units)
    angle_units = 'deg';
end

if nargin < 8 || isempty(dist_units)
    dist_units = 'm';
end

% Convert from LLA to ENU, response in meters
[e, n, u] = utils.lla2enu(lat, lon, alt, lat0, lon0, angle_units, dist_units);

[a, e, r] = utils.enu2aer(e, n, u, 'm');
