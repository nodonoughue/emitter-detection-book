function ell = loglikelihood(x_fdoa,v_fdoa,rho_dot,C,x_source,ref_idx)
% ell = loglikelihood(x_fdoa,v_fdoa,rho_dot,C,x_source,ref_idx)
%
% Computes the Log Likelihood for FDOA sensor measurement, given the 
% received measurement vector rho_dot, covariance matrix C, 
% and set of candidate source positions x_source.
%
% INPUTS:
%   x_fdoa      Sensor positions [m]
%   v_fdoa      Sensor velocities [m/s]
%   rho_dot     FDOA measurement vector
%   C           FDOA measurement error covariance matrix
%   x_source    Candidate source positions
%   ref_idx     Scalar index of reference sensor, or nDim x nPair
%               matrix of sensor pairings
%
% OUTPUTS:
%   ell         Log-likelihood evaluated at each position x_source.
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 6 || ~exist('ref_idx','var')
    ref_idx = [];
end

n_source_pos = size(x_source,2);
ell = zeros(n_source_pos,1);

% Resample covariance matrix
n_sensor = size(x_fdoa, 2);
[test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(ref_idx, n_sensor);
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
    r_dot = fdoa.measurement(x_fdoa,v_fdoa, x_i, ref_idx);
    
    % Evaluate the measurement error
    err = (rho_dot - r_dot);
    
    % Compute the scaled log likelihood
    if do_decomp
        ell(idx_source) = -err'/C_d*err;
    else
        ell(idx_source) = -err*C_inv*err;
    end
end