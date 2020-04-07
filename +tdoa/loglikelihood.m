function ell = loglikelihood(x_sensor,rho,C,x_source,ref_idx)
% ell = loglikelihood(x_sensor,rho,C,x_source,ref_idx)
%
% Computes the Log Likelihood for TDOA sensor measurement, given the 
% received measurement vector rho, covariance matrix C, 
% and set of candidate source positions x_source.
%
% INPUTS:
%   x_sensor    TDOA sensor positions 
%   rho         Measurement vector
%   C           measurement error covariance matrix
%   x_source    Candidate source position
%   ref_idx     Scalar index of reference sensor, or nDim x nPair
%               matrix of sensor pairings
%
% OUTPUTS:
%   ell         Log-likelihood evaluated at each position x_source.
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 5 || ~exist('ref_idx','var')
    ref_idx = [];
end

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
    x_i = x_source(:,idx_source);
    
    % Generate the ideal measurement matrix for this position
    r = tdoa.measurement(x_sensor, x_i, ref_idx);
    
    % Compute the measurement error
    err = (rho - r);

    % Evaluate the scaled log likelihood
    if do_decomp
        ell(idx_source) = -err'/C_d*err;
    else
        ell(idx_source) = -err*C_inv*err;
    end
end