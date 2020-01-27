function N = thermal_noise(bw,nf,t_ext)
% N = thermal_noise(bw,nf,t_ext)
%
% Compute the total noise power, given the receiver's noise bandwidth,
% noise figure, and external noise temperature.
%
% INPUTS:
%   bw      Receiver noise bandwidth [Hz]
%   nf      Receiver noise figure [dB] (DEFAULT = 0 dB)
%   t_ext   External noise temp [K] (DEFAULT = 0 K)
%
% OUTPUTS:
%   N       Total noise power [dBW]
%
% Nicholas O'Donoughue
% 1 July 2019

if nargin < 2 || isempty(nf)
    nf = 0;
end

if nargin < 3 || isempty(t_ext)
    t_ext = 0; % Standard temp
end

% Add the external noise temp to the reference temp (270 K)
T = utils.constants.T0 + t_ext;

% Boltzmann's Constant
k = 1.38e-23;

% Equation (D.6)
N = 10*log10(k.*T.*bw)+nf;