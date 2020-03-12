function Tm = moon_noise_temp()
% Tm = moon_noise_temp(f)
%
% Returns the noise temp (in Kelvin) for the moon at the specified
% frequency f (in Hertz). f can be a scalar, or N-dimensional matrix.
%
% INPUTS:
%   f           Carrier frequency [Hz]
%
% OUTPUTS:
%   Tm          Moon noise temp [K]
%
% Ref: Rec. ITU-R P.372-8
% Nicholas O'Donoughue
% 1 July 2019

% The moon noise temp is fairly constant across spectrum,
% with ~140 K during new moon phase and ~280 K during at full
% moon.  Using the arithmatic mean here as an approximate value.
Tm = (140 + 280)/2;