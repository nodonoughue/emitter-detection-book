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

n_dim = size(x_source,1);
n_source_pos = size(x_source,2);
ell = zeros(n_source_pos,1);

% Resample the covariance matrix
n_aoa = size(x_aoa,2);
n_tdoa = size(x_tdoa,2);
n_fdoa = size(x_fdoa,2);

% Determine if AOA measurements are 1D (azimuth) or 2D (az/el)
if n_aoa > 0
    assert(size(C,1) == n_aoa + n_tdoa + n_fdoa || size(C,1) == 2*n_aoa + n_tdoa + n_fdoa,'Unable to determine if AOA measurements are 1D or 2D');
    do2DAoA = size(C,1) == 2*n_aoa + n_tdoa + n_fdoa;
    if do2DAoA, n_aoa = 2*n_aoa; end
end

% Parse the TDOA and FDOA reference indices together
[tdoa_test_idx_vec, tdoa_ref_idx_vec] = utils.parseReferenceSensor(tdoa_ref_idx,n_tdoa);
[fdoa_test_idx_vec, fdoa_ref_idx_vec] = utils.parseReferenceSensor(fdoa_ref_idx,n_fdoa);
test_idx_vec = cat(2,1:n_aoa, n_aoa + tdoa_test_idx_vec, n_aoa + n_tdoa + fdoa_test_idx_vec);
ref_idx_vec = cat(2,nan(1,n_aoa), n_aoa + tdoa_ref_idx_vec, n_aoa + n_tdoa + fdoa_ref_idx_vec);

% For now, we assume the AOA is independent of TDOA/FDOA
C_tilde = utils.resampleCovMtx(C, test_idx_vec, ref_idx_vec);

% Ensure the covariance matrix is invertible
C_tilde = utils.ensureInvertible(C_tilde);

% Pre-compute covariance matrix inverses
do_decomp = ~verLessThan('MATLAB','9.3');
if do_decomp
    % Starging in R2017b, MATLAB released the DECOMPOSITION function,
    % which can decompose matrices for faster computation of left- and
    % right-division in for loops.
    C_d = decomposition(C_tilde,'chol');
else
    % If DECOMPOSITION is unavailable, let's precompute the pseudo-inverse.
    C_inv = pinv(C_tilde);
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