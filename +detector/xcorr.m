function detResult = xcorr(y1,y2,noise_var,M,prob_fa)
% detResult = xcorr(y1,y2,noise_var,M,prob_fa)
%
% Apply cross-correlation to determine whether a signal (y2) is present or
% absent in the provided received data vector (y1).
%
% Inputs:
%   y1          Data vector size (MxN)
%   y2          Desired signal (Mx1)
%   noise_var   Variance of the noise in y1
%   M           Number of samples in y1
%   prob_fa     Acceptable probability of false alarm
%
% Outputs:
%   detResult   Array of N binary detection results
%
% Nicholas O'Donoughue
% 1 July 2019

% Compute the sufficient statistic
sigma_0 = sqrt(M.*noise_var^2/2);
T = abs(sum(conj(y1).*y2,1)).^2./sigma_0.^2;

% Compute the threshold
if license('test','Statistics_Toolbox')
    eta = chi2inv(1-prob_fa,2);
else
    eta = utils.chi2inv(1-prob_fa,2);
end

% Compare T to eta
detResult = T > eta;

% In the rare event that T==eta, flip a weighted coin
coinFlipMask = T==eta;
detResult(coinFlipMask) = rand(sum(coinFlipMask),1)>(1-prob_fa);
