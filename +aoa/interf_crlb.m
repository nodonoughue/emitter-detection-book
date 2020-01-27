function crlb = interf_crlb(snr1,snr2,M,d_lam,psi_true)
% crlb = interf_crlb(snr1,snr2,M,d_lam,psi_true)
%
% Computes the lower bound on unbiased estimator error for an
% interferometer based direction of arrival receiver with multiple 
% amplitude samples taken (M samples)
%
% Inputs:
%
%   snr1        Signal-to-Noise ratio [dB] at receiver 1
%   snr2        Signal-to-Noise ratio [dB] at receiver 2
%   M           Number of samples
%   d_lam       Distance between receivers, divided by the signal
%               wavelength
%   psi_true    Angle of arrival [radians]
%
% Outputs:
%
%   crlb        Lower bound on the Mean Squared Error of an unbiased
%               estimation of psi (radians)
%
% Nicholas O'Donoughue
% 1 July 2019

% Compute the effective SNR
snr_lin1 = 10.^(snr1/10);
snr_lin2 = 10.^(snr2/10);
snr_eff = 1./(1./snr_lin1 + 1./snr_lin2);

crlb = (1./(2*M*snr_eff)).*(1./(2*pi*d_lam*cos(psi_true))).^2; % output in radians