function [x, y, z] = lla2ecef(lat, lon, alt, angle_units, dist_units)
% [x, y, z] = lla2ecef(lat, lon, alt, angle_units, dist_units)
%
% Convert a set of Lat, Lon, Alt coordinates to ECEF (x, y, z).  Lat,
% Lon, and Alt inputs must be broadcastable to a common size (all
% non-singleton dimensions must match).
%
% The optional inputs angle_units and dist_units are used to specify the
% units for lat/lon (either radians or degrees), and alt (any valid length
% unit).  See 'help unitsratio' for information on valid specifications.
%
% INPUTS:
%   lat             Vector or matrix of latitudes
%   lon             Vector or matrix of longitudes
%   alt             Vector or matrix of altitudes
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

if nargin < 4 || isempty(angle_units)
    angle_units = 'deg';
end

if nargin < 5 || isempty(dist_units)
    dist_units = 'm';
end

% Get scale factors
angle_mult = unitsratio('rad', angle_units);
dist_mult = unitsratio('m', dist_units);

% Lookup Earth Constants
semiMajorAxis_m = utils.constants.semimajor_axis_km*1e3;
ecc_sq = utils.constants.first_ecc_sq;

% Make sure Lat/Lon are in rad, alt in meters
lat_rad = lat * angle_mult;
lon_rad = lon * angle_mult;
alt_m = alt * dist_mult;

% Compute Sin and Cos of Lat/Lon Inputs
sin_lat = sin(lat_rad);
sin_lon = sin(lon_rad);
cos_lat = cos(lat_rad);
cos_lon = cos(lon_rad);

% Compute effective radius
eff_rad = semiMajorAxis_m ./ sqrt(1 - ecc_sq * sin_lat.^2);

% Compute ECEF Coordinates
x = (eff_rad + alt_m) .* cos_lat .* cos_lon;
y = (eff_rad + alt_m) .* cos_lat .* sin_lon;
z = ((1 - ecc_sq) * eff_rad + alt_m) .* sin_lat;