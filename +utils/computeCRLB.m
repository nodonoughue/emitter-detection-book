function crlb = computeCRLB(x,C,J,fim_prior,a_grad)
% crlb = computeCRLB(x,C,J)
%
% Generic form of the CRLB, valid for Gaussian distributed measurement
% models.
%
% Inputs:
%   x           (Ndim x M) array of source positions over which to 
%               calculate CRLB
%   C           Covariance matrix for measurements (N x N)
%   J           Jacobian matrix function handle; returns an Ndim x N matrix
%               representing the gradient of the log likelihood function
%               with respect to transmitter position. Accepts a single
%               input (Ndim x 1 transmitter position).
%   fim_prior   (Optional) Function handle that computes the Fisher 
%               Information Matrix for the statistical prior. If left
%               blank, or empty, then an uninformed prior is assumed.
%   a_grad      (Optional) Function handle that returns the gradient of the
%               equality constraints a(x) as an nDim x nConstraint 
%               matrix. If left blank, or empty, then no constraints are
%               applied.
%   
% Outputs:
%   crlb    Lower bound on the error covariance matrix for an unbiased
%           FDOA estimator (Ndim x Ndim)
%
% Nicholas O'Donoughue
% 16 May 2023

% Temporarily disable singular matrix warnings -- we know this happens
orig_state = warning('query','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:nearlySingularMatrix');

% Check inputs
do_prior = nargin >= 4 && ~isempty(fim_prior);
do_constraints = nargin >= 5 && ~isempty(a_grad);

% Parse Input Sizes
[n_dim, n_source] = size(x);

% Ensure the covariance matrix is invertible
C = utils.ensureInvertible(C);

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

% Initialize output variable
crlb = zeros([n_dim,n_dim,n_source]);

% Repeat CRLB for each of the n_source test positions
for idx =1:n_source
    this_x = x(:,idx);
    
    % Evaluate the Jacobian
    J_i = J(this_x);
    
    % Compute the Fisher Information Matrix
    if do_decomp
        F = J_i/C_d*J_i'; % Ndim x Ndim
    else
        F = J_i*C_inv*J_i';
    end
    
    % Evaluate the prior and add to the FIM
    if do_prior
        F = F + fim_prior(this_x);
    end

    if any(isnan(F(:))) || any(isinf(F(:)))
        % Problem is ill defined, Fisher Information Matrix cannot be
        % inverted
        crlb(:,:,idx) = NaN;
    else
        % Invert the FIM
        F_inv = pinv(F);
    
        if do_constraints
            % Evaluate the Gradient of the Constraint matrix
            A_i = a_grad(this_x);

            % Build the constraint effect matrix
            F_const_inv = F_inv * A_i / (A_i'*F_inv*A_i) * A_i'*F_inv;

            % CRLB is the inverse of the Fisher minus the constraint effect
            crlb(:,:,idx) = F_inv - F_const_inv;
        else
            crlb(:,:,idx) = F_inv;
        end
    end
end

% Restore warnings
warning(orig_state);