function Ta = atmospheric_noise(f,alt_start, el_angle_deg)
% Ta = atmospheric_noise(f,alt_start, el_angle_deg)
%
% Computes the noise temperature contribution from the reradaition of
% energy absorbed by the atmosphere in the direction of the antenna's 
% mainlobe.
%
% INPUTS:
%   f               Frequency [Hz]
%   alt_start       Altitude of receiver [m]
%   el_angle_deg    Elevation angle of receive mainbeam [degrees above
%                   local ground plane]
%
% OUTPUTS:
%   Ta               Atmospheric noise temperature [K]
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 2 || isempty(alt_start)
    alt_start = 0;
end

if nargin < 3 || isempty(el_angle_deg)
    el_angle_deg = 90;
end

% Assume integrated antenna gain is unity
alpha_a = 1;

% Compute zenith loss along main propagation path
zenith_angle_rad = (90-el_angle_deg)*pi/180;
L_dB = atm.calcZenithLoss(f,alt_start,zenith_angle_rad);
L = 10.^(L_dB/10);

% Compute average atmospheric temp
atmStruct = atm.standardAtmosphere(alt_start:100:100e3);
Tt=mean(atmStruct.T);
% Tt = utils.constants.T0;

% Equation D.12
Ta = alpha_a .* Tt .* (1-1./L);
