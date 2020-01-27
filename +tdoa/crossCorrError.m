function var = crossCorrError(snr,bw,pulseLen,bw_rms)
% var = crossCorrError(snr,bw,pulseLen,bw_rms)
%
% Computes the timing error for a Cross-Correlation time of arrival
% estimator, given the input signal's bandwidth, pulse length, and RMS
% bandwidth.
%
% Inputs:
%   snr         Signal-to-Noise ratio [dB]
%   bw          Bandwidth of input signal [Hz]
%   pulseLen    Length of input signal [s]
%   bw_rms      RMS Bandwidth of input signal [Hz]
%
% Outputs:
%
%   var         Timing error variance [s^2]
%
% Nicholas O'Donoughue
% 1 July 2019

% Convert input SNR to linear units
snr_lin = 10.^(snr/10);

% Compute the product of SNR, bandwidth, pulse length, and RMS bandwidth
A = bsxfun(@times,bsxfun(@times,snr_lin,bw),...
                  bsxfun(@times,pulseLen,bw_rms));

% Invert and apply 8*pi scale factor
var = 1./(8*pi*A);              