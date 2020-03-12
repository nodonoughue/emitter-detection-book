function Ts = sun_noise_temp(f)
% Ts = sun_noise_temp(f)
%
% Returns the noise temp (in Kelvin) for the sun at the specified
% frequency f (in Hertz). f can be a scalar, or N-dimensional matrix.
%
% Assumes a quiet sun, and represents a rough approximation from ITU
% documentation on radio noise.  Sun noise can be several orders of
% magnitude larger during solar disturbances.
%
% INPUTS:
%   f           Carrier frequency [Hz]
%
% OUTPUTS:
%   Ts          Sun noise temp [K]
%
% Ref: Rec. ITU-R P.372-14
% Nicholas O'Donoughue
% 1 July 2019

% Based on a visual reading on Figure 12 and the corresponding text
f_ghz = [.05,.2,1:10,20:10:100];
T_ref = [1e6,1e6,2e5,9e4,4.5e4,2.9e4,2e4,1.6e4,1.4e4,1.3e4,1.2e4,1e4,7e3,6.3e3,...
    6.2e3,6e3,6e3,6e3,6e3,6e3,6e3];

% Perform linear interpolation
Ts = interp1(f_ghz,T_ref,f/1e9,'linear',0);