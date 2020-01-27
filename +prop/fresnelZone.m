function fz = fresnelZone(f0,ht,hr)
% fz = fresnelZone(f0,ht,hr)
%
% Computes the Fresnel Zone for a given transmission,
% given by the equation:
%       FZ = 4*pi*h_t*h_r / lambda
%
% Inputs:
%   f0      Carrier frequency [Hz]
%   ht      Transmitter altitude [m]
%   hr      Receiver altitude [m]
%
% Outputs:
%   fz      Fresnel Zone range [m]
%
% Nicholas O'Donoughue
% 1 July 2019

% Convert the carrier frequency to wavelength [m]
lam = utils.constants.c./f0;

% Equation (B.3)
fz = 4*pi*ht.*hr./lam;