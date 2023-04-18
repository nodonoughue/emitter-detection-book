function Rh = radarHorizon(h1,h2)
% Rh = radarHorizon(h1,h2)
%
% Computes the radar horizon for transmitter and receiver at the given 
% distances above a smooth, round Earth.
%
% Leverages the (4/3) Earth radius approximation common for electromagnetic
% propagation.
%
% Inputs:
%   h1          Height of transmitter [m]
%   h2          Height of receiver [m]
%
% Outputs:
%   Rh          Radar horizon [m]
%
% Nicholas O'Donoughue
% 1 July 2019

R1 = sqrt(2*h1*utils.constant.Re+h1.^2);
R2 = sqrt(2*h2*utils.constant.Re+h2.^2);

Rh = bsxfun(@plus,R1,R2);