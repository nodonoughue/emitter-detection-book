function [lat, lon, alt] = aer2lla(a, e, r, lat0, lon0, alt0, angle_units, dist_units)
% [lat, lon, alt] = aer2lla(a, e, r, angle_units, dist_units)
%
% Convert spherical AER coordinates to LLA, where azimuth is
% in degrees East from North (as opposed to the typical degrees +y from +x)
% and elevation is degrees above the horizontal.
%
% INPUTS:
%   a               Azimuth (deg E from N)
%   e               Elevation (deg above horizontal)
%   r               Range (meters)
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

% Convert to ECEF
[x, y, z] = utils.aer2ecef(a, e, r, lat0, lon0, alt0, angle_units, dist_units);

% Conver to LLA
[lat, lon, alt] = utils.ecef2lla(x, y, z, 'm');