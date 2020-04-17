function xi = xcorrMinSNR(PFA,PD,Tcorr,Tp,Bn,Bs)
% xi = xcorrMinSNR(PFA,PD,Tcorr,Tp,Bn,Bs)
%
% Compute the required SNR to achieve the desired probability of detection,
% given the maximum acceptable probability of false alarm, and the number
% of complex samples M.
%
% The returned SNR is the ratio of signal power to complex noise power.
%
% Inputs:
%
%   PFA             Probability of False Alarm [0-1]
%   PD              Probability of Detection [0-1]
%   Tcorr           Correlation time [sec]
%   Tp              Pulse Duration [sec]
%   Bn              Noise bandwidth [Hz]
%   Bs              Signal Bandwidth [Hz]
%
% Outputs
%   xi              Signal-to-Noise ratio [dB]
%
% Nicholas O'Donoughue
% 1 July 2019

% Make sure the signal bandwidth and time are observable
Bs = min(Bs,Bn);
Tp = min(Tp,Tcorr);
M = fix(Tcorr*Bn);

% Find the min SNR after cross-correlation processing
xi_min_out = detector.squareLawMinSNR(PFA,PD,M);

% Invert the SNR Gain equation
%   xi_out = 2*T*B*xi_in^2 / 1+(2*xi_in)
xi_out_lin = 10.^(xi_min_out/10);
xi_in_lin = (xi_out_lin + sqrt(xi_out_lin*(xi_out_lin+Bs*Tcorr)))/(Tp*Bs);
xi = 10*log10(xi_in_lin);


% Method 2
% if license('test','Statistics_Toolbox')
%     eta = chi2inv(1-PFA,2); % 2 DOF
%      pd_check = 1-ncx2cdf(eta/(1+2*xi_in_lin),2,xi_out_lin);
% else
%     eta = utils.chi2inv(1-PFA,2);
%      pd_check = 1-utils.ncx2cdf(eta/(1+2*xi_in_lin),2,xi_out_lin);
% end