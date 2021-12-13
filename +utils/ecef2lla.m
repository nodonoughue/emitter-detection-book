function [lat, lon, alt] = ecef2lla(x, y, z, dist_units)
% [lat, lon, alt] = ecef2lla(x, y, z, dist_units)
%
% Conversion from cartesian ECEF to geodetic LLA coordinates, using a
% direct computation methods.  Note that iterative solutions also exist,
% and can be made to be more accurate (by iterating until a desired
% accuracy is achieved).  For the purposes of this package, we have
% implemented the faster direct computation method, understanding that
% its accuracy may be limited in certain edge cases.
%
%
% INPUTS:
%   x               Vector or matrix of ECEF x-coordinates
%   y               Vector or matrix of ECEF y-coordinates
%   z               Vector or matrix of ECEF z-coordinates
%   dist_units      Units for ECEF coordinates [Default='m']
%
% OUTPUTS:
%   lat             Latitude (degrees North)
%   lon             Longitude (degrees East)
%   alt             Altitude (meters)
%
% Source: https://microem.ru/files/2012/08/GPS.G1-X-00006.pdf
%
% Nicholas O'Donoughue
% 9 June 2021

if nargin < 4 || isempty(dist_units)
    dist_units = 'm';
end

% Compute the longitude as the arctan of y/x
lon_rad = atan2(y, x);

% Read in constants
a = utils.constants.semimajor_axis_km * 1e3;
b = utils.constants.semiminor_axis_km * 1e3;
e1_sq = utils.constants.first_ecc_sq;
e2_sq = utils.constants.second_ecc_sq;

% Compute Auxiliary Values
p = sqrt(x.^2 + y.^2);
th= atan2(z*a, p*b);

% Compute Latitude estimate
lat_rad = atan2(z + e2_sq * b * sin(th).^3, p - e1_sq*a*cos(th).^3);

% Compute Height
N = a ./ sqrt(1-e1_sq*sin(lat_rad).^2);
alt_m = p./cos(lat_rad) - N;

% Format Outputs
lat = lat_rad * 180/pi;
lon = lon_rad * 180/pi;
alt = alt_m;