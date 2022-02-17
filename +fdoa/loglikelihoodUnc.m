function ell = loglikelihoodUnc(x_fdoa, v_fdoa, zeta, C, C_beta, theta, fdoa_ref_idx)
% function ell = loglikelihoodUnc(x_fdoa, v_fdoa, zeta, C, C_beta, ...
%                                 theta, fdoa_ref_idx)
%
% Computes the Log Likelihood for FDOA, given the received measurement 
% vector zeta, covariance matrix C, and parameter vector 
% theta = [x; alpha; beta], where alpha is the set of measurement biases, 
% and beta is the set of sensor positions and covariances.
%
% INPUTS:
%   x_fdoa      nDim x nFDOA vector of reported FDOA sensor positions
%   v_fdoa      nDim x nFDOA vector of reported FDOA sensor velocities
%   zeta        Combined AOA/TDOA/FDOA measurement vector
%   C           Combined AOA/TDOA/FDOA measurement covariance matrix
%   C_beta      nDim x (nAOA + nTDOA + 2*nFDOA) sensor pos/vel error
%               covariance matrix
%   theta       Candidate parameter vectors theta=[x;alpha;beta]
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA measurements
%
% OUTPUTS:
%   ell         Log-likelihood evaluated at each position x_source.
%
% Nicholas O'Donoughue
% 22 Feb 2022

% Parse inputs
if nargin < 7 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

%% Parse Dimensions
n_sets = size(theta,2);
ell = zeros(n_sets,1);

n_dim = size(x_fdoa,1);
n_fdoa = size(x_fdoa,2);

% Determine how many TDOA/FDOA measurements there are
[fdoa_test_idx_vec, fdoa_ref_idx_vec] = utils.parseReferenceSensor(fdoa_ref_idx,n_fdoa);

% Assemble mean value for beta; given by the reported sensor positions and
% velocities
beta0 = [x_fdoa(:); v_fdoa(:)];

%% Compute reference indices for the elements of theta
% theta vector contains x, alpha, and beta.  Let's define the
% indices
x_ind = 1:n_dim;
alpha_f_ind = n_dim + (1:n_fdoa);
beta_fx_ind = alpha_f_ind(end) + (1:n_dim*n_fdoa);
beta_fv_ind = beta_fx_ind(end)+ (1:n_dim*n_fdoa);

%% Resample the covariance matrix

% For now, we assume the AOA is independent of TDOA/FDOA
C_tilde = utils.resampleCovMtx(C, fdoa_test_idx_vec, fdoa_ref_idx_vec);

% Ensure the covariance matrix is invertible
C_tilde = utils.ensureInvertible(C_tilde);

% Pre-compute covariance matrix inverses
do_decomp = ~verLessThan('MATLAB','9.3');
if do_decomp
    % Starging in R2017b, MATLAB released the DECOMPOSITION function,
    % which can decompose matrices for faster computation of left- and
    % right-division in for loops.
    C_d = decomposition(C_tilde,'chol');
    C_b = decomposition(C_beta,'chol');
else
    % If DECOMPOSITION is unavailable, let's precompute the pseudo-inverse.
    C_inv = pinv(C_tilde);
    C_b_inv = pinv(C_beta);
end



for idx_source = 1:n_sets
    th_i = theta(:,idx_source);
    x_i = th_i(x_ind);
    alpha_fdoa = th_i(alpha_f_ind);
    x_fdoa_i = reshape(th_i(beta_fx_ind),n_dim,n_fdoa);
    v_fdoa_i = reshape(th_i(beta_fv_ind),n_dim,n_fdoa);
    beta_i = [x_fdoa_i(:); v_fdoa_i(:)];

    % Generate the ideal measurement matrix for this position
    z = fdoa.measurement(x_fdoa,v_fdoa, x_i, fdoa_ref_idx, alpha_fdoa);
    
    % Evaluate the measurement error
    err = (z - zeta);
    err_b = (beta0 - beta_i);

    % Evaluate the 
    % Compute the scaled log likelihood
    if do_decomp
        ell(idx_source) = -err'/C_d*err - err_b'/C_b*err_b;
    else
        ell(idx_source) = -err'*C_inv*err - err_b'*C_b_inv*err_b;
    end
end