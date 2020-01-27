function ell = loglikelihood(x_aoa, psi,C,x_source)
% function ell = loglikelihood(x_aoa,psi,C,x_source)
%
% Computes the Log Likelihood for AOA sensor measurement, given the 
% received measurement vector psi, covariance matrix C, 
% and set of candidate source positions x_source.
%
% INPUTS:
%   x_aoa       nDim x nAOA vector of AOA sensor positions
%   psi         Measurement vector
%   C           Measurement error covariance matrix
%   x_source    Candidate source positions
%
% OUTPUTS:
%   ell         Log-likelihood evaluated at each position x_source.
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
n_source_pos = size(x_source,2);
ell = zeros(n_source_pos,1);

% Preprocess covariance matrix
C_d = decomposition(C);

for idx_source = 1:n_source_pos
    % Generate the ideal measurement matrix for this position
    p = triang.measurement(x_aoa, x_source(:,idx_source));
    
    % Compute the log-likelihood
    err = (psi - p);
    
    ell(idx_source) = -err'/C_d*err;
end