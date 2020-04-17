function detResult = squareLaw(z,noise_var,prob_fa)
% detResult = squareLaw(z,noise_var,prob_fa)
%
% Compute detection via the square law and return binary detection events
%
% Inputs:
%   z               Input signal, (MxN) for M samples per detection event,
%                   and N separate test
%   noise_var       Noise variance on input signal
%   prob_fa         Acceptable probability of false alarm
%
% Outputs:
%   detResult       Array of N binary detection results
%
% Nicholas O'Donoughue
% 1 July 2019

% Compute the sufficient statistic
T = sum(abs(z).^2,1)/noise_var;

% Compute the threshold
if license('test','Statistics_Toolbox')
    eta = chi2inv(1-prob_fa,2*size(z,1));
else
    eta = utils.chi2inv(1-prob_fa,2*size(z,1));
end

% Compare T to eta
detResult = T > eta;

% In the rare event that T==eta, flip a weighted coin
coinFlipMask = T==eta;
detResult(coinFlipMask) = rand(sum(coinFlipMask),1)>(1-prob_fa);
