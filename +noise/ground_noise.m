function Tg = ground_noise(G_gnd_dB, E, angular_area)
% Tg = ground_noise(G_gnd_dB, E, angular_area)
%
% Compute the combined noise temperature from ground effects; predominantly 
% caused by reradiation of thermal energy from the sun.
%
% INPUTS:
%   G_gnd_dB        Average antenna gain in direction of the ground [dBi]
%                   (DEFAULT = -5 dBi)
%   E               Emissivity of ground (Default = 1)
%   angular_area    Area (in steradians) of ground as visible from antenna,
%                   (DEFAULT = pi)
% 
% OUTPUTS:
%   Tg      Ground noise temperature [K]
%
% Nicholas O'Donoughue
% 1 July 2019

if nargin < 1 || isempty(G_gnd_dB)
    G_gnd_dB = -5;
end

if nargin < 2 || isempty(E)
    % Assume Emissivity is 1
    E = 1;
end

if nargin < 3 || isempty(angular_area)
    % Assume pi steradians effective angular
    % area of ground as seen by antenna
    angular_area = pi;
end

% Convert average ground antenna gain to linear units
G_gnd = 10.^(G_gnd_dB/10);

% Assume ground temp is 290 K
Tt = 290;

% Compute ground noise temp according to (D.13)
Tg = angular_area .* G_gnd .* E .* Tt / (4*pi);
