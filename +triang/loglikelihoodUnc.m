function ell = loglikelihoodUnc(x_aoa, zeta, C, C_beta, theta)
% function ell = loglikelihoodUnc(x_aoa, zeta, C, C_beta, theta)
%
% Computes the Log Likelihood for AOA, given the received measurement 
% vector zeta, covariance matrix C, and parameter vector 
% theta = [x; alpha; beta], where alpha is the set of measurement biases, 
% and beta is the set of sensor positions and covariances.
%
% INPUTS:
%   x_aoa       nDim x nAOA vector of reported AOA sensor positions
%   zeta        Combined AOA/TDOA/FDOA measurement vector
%   C           Combined AOA/TDOA/FDOA measurement covariance matrix
%   C_beta      nDim x (nAOA + nTDOA + 2*nFDOA) sensor pos/vel error
%               covariance matrix
%   theta       Candidate parameter vectors theta=[x;alpha;beta]
%
% OUTPUTS:
%   ell         Log-likelihood evaluated at each position x_source.
%
% Nicholas O'Donoughue
% 22 Feb 2022

%% Parse Dimensions
n_sets = size(theta,2);
ell = zeros(n_sets,1);

n_dim = size(x_aoa,1);
n_aoa = size(x_aoa,2);

% Determine if AOA measurements are 1D (azimuth) or 2D (az/el)
assert(size(C,1) == n_aoa || size(C,1) == 2*n_aoa,'Unable to determine if AOA measurements are 1D or 2D');
do2DAoA = size(C,1) == 2*n_aoa + n_tdoa + n_fdoa;
if do2DAoA
    m_aoa = 2*n_aoa;
else
    m_aoa = n_aoa;
end

% Assemble mean value for beta; given by the reported sensor positions and
% velocities
beta0 = x_aoa(:);

%% Compute reference indices for the elements of theta
% theta vector contains x, alpha, and beta.  Let's define the
% indices
x_ind = 1:n_dim;
alpha_ind = x_ind(end) + (1:m_aoa);
beta_ind = alpha_ind(end) + (1:n_dim*n_aoa);

%% Resample the covariance matrix

% Ensure the covariance matrix is invertible
C = utils.ensureInvertible(C);

% Pre-compute covariance matrix inverses
do_decomp = ~verLessThan('MATLAB','9.3');
if do_decomp
    % Starging in R2017b, MATLAB released the DECOMPOSITION function,
    % which can decompose matrices for faster computation of left- and
    % right-division in for loops.
    C_d = decomposition(C,'chol');
    C_b = decomposition(C_beta,'chol');
else
    % If DECOMPOSITION is unavailable, let's precompute the pseudo-inverse.
    C_inv = pinv(C);
    C_b_inv = pinv(C_beta);
end



for idx_source = 1:n_sets
    th_i = theta(:,idx_source);
    x_i = th_i(x_ind);
    alpha_aoa = th_i(alpha_ind);
    x_aoa_i = reshape(th_i(beta_ind),n_dim,n_aoa);
    beta_i = x_aoa_i(:);

    % Generate the ideal measurement matrix for this position
    z = triang.measurement(x_aoa, x_i, do2DAoA, alpha_aoa);
    
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