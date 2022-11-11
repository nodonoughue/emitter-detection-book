function ell = loglikelihoodUnc(x_aoa, x_tdoa, x_fdoa, v_fdoa, zeta, C, C_beta, theta, tdoa_ref_idx, fdoa_ref_idx)
% function ell = loglikelihoodUnc(x_aoa, x_tdoa, x_fdoa, v_fdoa, zeta, ...
%                           C, C_beta, theta, tdoa_ref_idx, fdoa_ref_idx)
%
% Computes the Log Likelihood for Hybrid sensor measurement (AOA, TDOA, and
% FDOA), given the received measurement vector zeta, covariance matrix C, 
% and parameter vector theta = [x; alpha; beta], where alpha is the set of
% measurement biases, and beta is the set of sensor positions and
% covariances.
%
% INPUTS:
%   x_aoa       nDim x nAOA vector of reported AOA sensor positions
%   x_tdoa      nDim x nTDOA vector of reported TDOA sensor positions
%   x_fdoa      nDim x nFDOA vector of reported FDOA sensor positions
%   v_fdoa      nDim x nFDOA vector of reported FDOA sensor velocities
%   zeta        Combined AOA/TDOA/FDOA measurement vector
%   C           Combined AOA/TDOA/FDOA measurement covariance matrix
%   C_beta      nDim x (nAOA + nTDOA + 2*nFDOA) sensor pos/vel error
%               covariance matrix
%   theta       Candidate parameter vectors theta=[x;alpha;beta]
%   tdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for TDOA measurements
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA measurements
%
% OUTPUTS:
%   ell         Log-likelihood evaluated at each position x_source.
%
% Nicholas O'Donoughue
% 22 Feb 2022

% Parse inputs
if nargin < 8 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

if nargin < 9 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

%% Parse Dimensions
n_sets = size(theta,2);
ell = zeros(n_sets,1);

n_dim = size(x_aoa,1);
n_aoa = size(x_aoa,2);
n_tdoa = size(x_tdoa,2);
n_fdoa = size(x_fdoa,2);

% Determine if AOA measurements are 1D (azimuth) or 2D (az/el)
assert(size(C,1) == n_aoa + n_tdoa + n_fdoa || size(C,1) == 2*n_aoa + n_tdoa + n_fdoa,'Unable to determine if AOA measurements are 1D or 2D');
do2DAoA = size(C,1) == 2*n_aoa + n_tdoa + n_fdoa;
if do2DAoA
    m_aoa = 2*n_aoa;
else
    m_aoa = n_aoa;
end

% Determine how many TDOA/FDOA measurements there are
[tdoa_test_idx_vec, tdoa_ref_idx_vec] = utils.parseReferenceSensor(tdoa_ref_idx,n_tdoa);
[fdoa_test_idx_vec, fdoa_ref_idx_vec] = utils.parseReferenceSensor(fdoa_ref_idx,n_fdoa);

% Assemble mean value for beta; given by the reported sensor positions and
% velocities
beta0 = [x_aoa(:); x_tdoa(:); x_fdoa(:); v_fdoa(:)];

%% Compute reference indices for the elements of theta
% theta vector contains x, alpha, and beta.  Let's define the
% indices
x_ind = 1:n_dim;
alpha_a_ind = x_ind(end) + (1:m_aoa);
alpha_t_ind = alpha_a_ind(end) + (1:n_tdoa);
alpha_f_ind = alpha_t_ind(end) + (1:n_fdoa);
beta_a_ind = alpha_f_ind(end) + (1:n_dim*n_aoa);
beta_t_ind = beta_a_ind(end) + (1:n_dim*n_tdoa);
beta_fx_ind = beta_t_ind(end) + (1:n_dim*n_fdoa);
beta_fv_ind = beta_fx_ind(end)+ (1:n_dim*n_fdoa);

%% Resample the covariance matrix

% Parse the TDOA and FDOA reference indices together
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
    C_b = decomposition(C_beta,'chol');
else
    % If DECOMPOSITION is unavailable, let's precompute the pseudo-inverse.
    C_inv = pinv(C_tilde);
    C_b_inv = pinv(C_beta);
end



for idx_source = 1:n_sets
    th_i = theta(:,idx_source);
    x_i = th_i(x_ind);
    alpha_aoa = th_i(alpha_a_ind);
    alpha_tdoa = th_i(alpha_t_ind);
    alpha_fdoa = th_i(alpha_f_ind);
    x_aoa_i = reshape(th_i(beta_a_ind),n_dim,n_aoa);
    x_tdoa_i = reshape(th_i(beta_t_ind),n_dim,n_tdoa);
    x_fdoa_i = reshape(th_i(beta_fx_ind),n_dim,n_fdoa);
    v_fdoa_i = reshape(th_i(beta_fv_ind),n_dim,n_fdoa);
    beta_i = [x_aoa_i(:); x_tdoa_i(:); x_fdoa_i(:); v_fdoa_i(:)];

    % Generate the ideal measurement matrix for this position
    z = hybrid.measurement(x_aoa,x_tdoa,x_fdoa,v_fdoa, x_i, tdoa_ref_idx, fdoa_ref_idx, do2DAoA, alpha_aoa, alpha_tdoa, alpha_fdoa);
    
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