function r_eff = effRadiusEarth(lat, is_deg)
% r_eff = effRadiusEarth(lat, is_deg)
%
% Computes the effective radius of the Earth, at a given latitud (lat).
%
% The latitude is geodetic (not geocentric), the result is computed 
% according to equation 4.13:
%     r_eff = a / sqrt(1-e_1^2 sin^2(lat))
%
% where a = semi-major axis, in meters, and e_1^2 is the square of the
% first eccentricity of the Earth's ellipsoid.
%
% INPUTS:
%   lat         Geodetic latitude (deg/rad)
%   is_deg      [Optional] flag indicating whether input is in degrees or
%               radians (default=True -- degrees)
%
% OUTPUTS:
%   r_eff       Effective radius of the Earth at the specified Latitude
%
% Nicholas O'Donoughue
% 27 August 2021

% Load Constants
a = utils.constants.semimajor_axis_km * 1e3;
e1sq = utils.constants.first_ecc_sq;

% Parse Inputs
if nargin < 2 || isempty(is_deg)
    is_deg = true;
end

if is_deg
    lat_rad = lat*pi/180;
else
    lat_rad = lat;
end

% Equation 4.13
r_eff = a ./ sqrt(1 - e1sq*sin(lat).^2);