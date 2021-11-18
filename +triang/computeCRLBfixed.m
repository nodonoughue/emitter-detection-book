function crlb = computeCRLBfixed(x_aoa,x0,C,a_grad)
% crlb = computeCRLBfixed(x_aoa,xs,C,a_grad)
%
% Computes the CRLB on position accuracy for source at location xs and
% sensors at locations in x_aoa (Ndim x N).  C is an NxN matrix of TOA
% covariances at each of the N sensors.
%
% Employs the constrained CRLB, according to equation 5.15.
%
% Inputs:
%   x_aoa           (Ndim x N) array of AOA sensor positions
%   x0              (Ndim x M) array of source positions over which to 
%                   calculate CRLB
%   C               AOA covariance matrix (radians^2)
%   a_grad          Function handle that returns the gradient of the
%                   equality constraints a(x) as an nDim x nConstraint 
%                   matrix
%
% Outputs:
%   crlb    Lower bound on the error covariance matrix for an unbiased
%           AOA estimator (Ndim x Ndim)
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
n_dim = size(x_aoa,1);
n_source = size(x0,2);

% Set up Jacobian function
J = @(x) triang.jacobian(x_aoa,x);

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

warning('off','MATLAB:nearlySingularMatrix');

% Repeat CRLB for each of the n_source test positions
for idx =1:n_source
    this_x = x0(:,idx);
    
    % Compute Jacobian matrix
    J_i = J(this_x);
    
    % Evaluate the Gradient of the Constraint matrix
    A_i = a_grad(this_x);

    % Compute Fisher Information Matrix
    if do_decomp
        F = J_i/C_d*J_i'; % Ndim x Ndim
    else
        F = J_i*C_inv*J_i';
    end
    
    if any(isnan(F(:))) || any(isinf(F(:)))
        % Problem is ill defined, Fisher Information Matrix cannot be
        % inverted
        crlb(:,:,idx) = NaN;
        continue;
    else    
        % Invert the Fisher Information Matrix to compute the CRLB
        F_inv = pinv(F);

        % Build the constraint effect matrix
        F_const_inv = F_inv * A_i / (A_i'*F_inv*A_i) * A_i'*F_inv;

        % CRLB is the inverse of the Fisher minus the constraint effect
        crlb(:,:,idx) = F_inv - F_const_inv;
    end
end

warning('on','MATLAB:nearlySingularMatrix');