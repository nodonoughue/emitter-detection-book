function sigma = freqEstCRLB(sample_time,num_samples,snr_db)
% sigma = freqEstCRLB(time_s,bw_hz,snr_db)
%
% Compute the CRLB for the frequency difference estimate from a pair of
% sensors, given the time duration of the sampled signals, receiver
% bandwidth, and average SNR.
%
% INPUTS:
%   sample_time Received signal duration [s]
%   num_samples Number of receiver samples
%   snr_db      SNR [dB]
%
% OUTPUTS:
%   sigma       Frequency difference estimate error standard deviation [Hz]
%
% Nicholas O'Donoughue
% 1 July 2019

% Convert SNR to linear units
snr_lin = 10.^(snr_db/10);

% Compute the CRLB of the center frequency variance
sigma = sqrt(3./(pi.^2*sample_time.^2.*num_samples.*(num_samples.^2-1).*snr_lin));
