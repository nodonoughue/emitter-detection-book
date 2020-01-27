function err = peakDetectionError(snr)
% err = peakDetectionError(snr)
%
% Computes the error in time of arrival estimation for a peak detection 
% algorithm, based on input SNR.
%
% INPUTS
%   snr         Signal-to-Noise Ratio [dB]
%
% OUTPUTS
%
%   err         Expected error variance
%
% Nicholas O'Donoughue
% 1 July 2019

% Convert SNR to linear units
snr_lin = 10.^(snr/10);

% Compute Error
err = 1./(2*snr_lin);