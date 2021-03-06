function ell = loglikelihood(x_aoa,x_tdoa,x_fdoa,v_fdoa,zeta,C,x_source,tdoa_ref_idx,fdoa_ref_idx)
% function ell = loglikelihood(x_aoa,x_tdoa,x_fdoa,zeta,C,x_source,...
%                                               tdoa_ref_idx,fdoa_ref_idx)
%
% Computes the Log Likelihood for Hybrid sensor measurement (AOA, TDOA, and
% FDOA), given the received measurement vector zeta, covariance matrix C, 
% and set of candidate source positions x_source.
%
% INPUTS:
%   x_aoa       nDim x nAOA vector of AOA sensor positions
%   x_tdoa      nDim x nTDOA vector of TDOA sensor positions
%   x_fdoa      nDim x nFDOA vector of FDOA sensor positions
%   v_fdoa      nDim x nFDOA vector of FDOA sensor velocities
%   zeta        Combined AOA/TDOA/FDOA measurement vector
%   C           Combined AOA/TDOA/FDOA measurement covariance matrix
%   x_source    Candidate source positions
%   tdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for TDOA measurements
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA measurements
%
% OUTPUTS:
%   ell         Log-likelihood evaluated at each position x_source.
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 8 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

if nargin < 9 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
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
    z = hybrid.measurement(x_aoa,x_tdoa,x_fdoa,v_fdoa, x_i, tdoa_ref_idx, fdoa_ref_idx);
    
    % Evaluate the measurement error
    err = (z - zeta);

    % Compute the scaled log likelihood
    if do_decomp
        ell(idx_source) = -err'/C_d*err;
    else
        ell(idx_source) = -err'*C_inv*err;
    end
end