function ell = loglikelihood(x_sensor, psi,C,x_source)
% function ell = loglikelihood(x_sensor,psi,C,x_source)
%
% Computes the Log Likelihood for AOA sensor measurement, given the 
% received measurement vector psi, covariance matrix C, 
% and set of candidate source positions x_source.
%
% INPUTS:
%   x_sensor    nDim x nAOA vector of AOA sensor positions
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

% Pre-compute covariance matrix inverses
do_decomp = ~verLessThan('MATLAB','9.3');
if do_decomp
    % Starging in R2017b, MATLAB released the DECOMPOSITION function,
    % which can decompose matrices for faster computation of left- and
    % right-division in for loops.
    C_d = decomposition(C,'chol');
else
    % If DECOMPOSITION is unavailable, let's precompute the pseudo-inverse.
    C_inv = pinv(C);
end
for idx_source = 1:n_source_pos
    % Generate the ideal measurement matrix for this position
    p = triang.measurement(x_sensor, x_source(:,idx_source));
    
    % Compute the measurement error
    err = (psi - p);

    % Evaluate the scaled log likelihood
    if do_decomp
        ell(idx_source) = -err'/C_d*err;
    else
        ell(idx_source) = -err*C_inv*err;
    end
end