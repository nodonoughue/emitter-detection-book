function crlb = watson_watt_crlb(snr,M)
% crlb = watson_watt_crlb(snr,M)
%
% Compute the lower bound on unbiased estimation error for a Watson-Watt
% based angle of arrival receiver.
%
% Inputs:
%
%   snr         Signal-to-Noise ratio [dB]
%   M           The number of samples taken
%
% Outputs:
%
%   crlb        Lower bound on the Mean Squared Error of an unbiased
%               estimation of psi (radians)
%
% Nicholas O'Donoughue
% 1 July 2019


snr_lin = 10.^(snr/10);
crlb = 1./(M.*snr_lin);