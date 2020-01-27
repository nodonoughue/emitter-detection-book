function crlb = directional_crlb(snr,M,g,g_dot,psi_samples,psi_true)
% crlb = directional_crlb(snr,M,g,g_dot,psi_samples,psi_true)
%
% Computes the CRLB for a directional antenna with amplitude measurements
% taken at a series of angles.  Supports M measurements from each of N
% different angles.
%
% If there are multiple true angles of arrival provided (psi_true), then
% the CRLB is computed independently for each one.
%
% Inputs:
%
%   snr         Signal-to-Noise ratio [dB]
%   M           Number of samples for each antenna position
%   g           Function handle to g(psi)
%   g_dot       Function handle to g_dot(psi)
%   psi_samples The sampled steering angles (radians)
%   psi_true    The true angle of arrival (radians)
%
% Outputs:
%
%   crlb        Lower bound on the Mean Squared Error of an unbiased
%               estimation of psi (radians)
%
% Nicholas O'Donoughue
% 1 July 2019

% Convert SNR from dB to linear units
snr_lin = 10.^(snr/10);

% Evaluate the antenna pattern and antenna gradient at each of the steering
% angles sampled.
g_vec = arrayfun(@(x) g(x),psi_samples(:)-psi_true);
g_dot_vec = arrayfun(@(x) g_dot(x),psi_samples(:)-psi_true);

% Pre-compute steering vector inner products
g_g = sum(g_vec.*g_vec,1);
g_dot_g = sum(g_vec.*g_dot_vec,1);
g_dot_g_dot = sum(g_dot_vec.*g_dot_vec,1);

% Compute CRLB for each true angle theta
J = 2*M*snr_lin*(g_dot_g_dot-(g_dot_g).^2./g_g);
crlb = 1./J; % 1 x numel(th)