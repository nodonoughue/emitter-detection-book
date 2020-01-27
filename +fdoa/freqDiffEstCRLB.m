function sigma = freqDiffEstCRLB(time_s,bw_hz,snr_db)
% sigma = freqDiffEstCRLB(time_s,bw_hz,snr_db)
%
% Compute the CRLB for the frequency difference estimate from a pair of
% sensors, given the time duration of the sampled signals, receiver
% bandwidth, and average SNR.
%
% INPUTS:
%   time_s      Received signal duration [s]
%   bw_hz       Received signal bandwidth [Hz]
%   snr_db      Average SNR [dB]
%
% OUTPUTS:
%   sigma       Frequency difference estimate error standard deviation [Hz]
%
% Nicholas O'Donoughue
% 1 July 2019

% Convert SNR to linear units
snr_lin = 10.^(snr_db/10);

% Apply the CRLB equations
sigma = sqrt(3 ./ (4.*pi.^2*time_s.^3.*bw_hz.*snr_lin));