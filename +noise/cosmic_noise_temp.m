function Tc = cosmic_noise_temp(f,rx_alt,alpha_c,G_sun_dB,G_moon_dB)
% Tc = cosmic_noise_temp(f,rx_alt,alpha_c,G_sun_dB,G_moon_dB)
%
% Computes the combined cosmic noise temperature, including contributions
% from the sun, the moon, and the galactic background.  Includes
% approximate effect of atmospheric loss (sun and moon are treated as
% as coming from zenith; rather than their true angles.
%
% INPUTS:
%   f           Carrier frequency [Hz]
%   rx_alt      Receiver altitude [m]
%   alpha_c     Fraction of the antenna's receive pattern that is above
%               the horizon [0-1]   
%   G_sun_dB    Antenna gain directed at the sun [dBi]
%   G_moon_dB   Antenna gain directed at the moon [dBi]
%
% OUTPUTS:
%   Tc          Combined cosmic noise temperature [K]
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse Inputs
if nargin < 2 || isempty(rx_alt)
    rx_alt = 0;
end

if nargin < 3 || isempty(alpha_c)
    alpha_c = .95;
end

if nargin < 4 || isempty(G_sun_dB)
    G_sun_dB = -Inf;
end

if nargin < 5 || isempty(G_moon_dB)
    G_moon_dB = -Inf;
end

%% Compute Raw Noise Temp
T_100 = 3050; % Geometric mean of 100 MHz noise spectrum samples

Tc0 = T_100 * (100e6./f).^(2.5) + 2.7;
Ts0 = noise.sun_noise_temp(f);
Tm0 = noise.moon_noise_temp();

% Above 2 GHz, the only contribution is from cosmic background radiation
% (2.7 K), which is essentially negligible.
freq_mask = f < 2e9;
Tc0(~freq_mask) = 2.7;


%% Apply Antenna Patterns
G_sun = 10.^(G_sun_dB/10);
G_moon = 10.^(G_moon_dB/10);

T0 = Tc0 * alpha_c + ...
     Ts0 * 4.75e-6 * G_sun + ...
     Tm0 * 4.75e-6 * G_moon;

%% Apply Atmospheric Losses

% Get zenith loss
LA_dB = reshape(atm.calcZenithLoss(f,rx_alt,pi/4),size(f));
LA = 10.^(LA_dB/10);

Tc = T0./LA;