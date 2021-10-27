function [a, e, r] = ecef2aer(x, y, z, lat0, lon0, alt0, angle_units, dist_units)
% [a, e, r] = ecef2aer(x, y, z, lat0, lon0, angle_units, dist_units)
%
% Convert a set of ECEF coordinates to AER (azimuth, elevation, range),
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

% Convert to ENU
[e, n, u] = utils.ecef2enu(x, y, z, lat0, lon0, alt0, angle_units, dist_units);

% Convert ENU to AER
[a, e, r] = utils.enu2aer(e, n, u, 'm');