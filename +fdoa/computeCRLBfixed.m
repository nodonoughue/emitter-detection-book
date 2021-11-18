function crlb = computeCRLBfixed(x_fdoa,v_fdoa,xs,C,a_grad,ref_idx)
% crlb = computeCRLBfixed(x_fdoa,v_fdoa,xs,C,a_grad,ref_idx)
%
% Computes the CRLB on position accuracy for source at location xs and
% sensors at locations in x_fdoa (Ndim x N) with velocity v_fdoa.  
% C is an Nx1 vector of FOA variances at each of the N sensors, and ref_idx
% defines the reference sensor(s) used for FDOA.
%
% Employs the constrained CRLB, according to equation 5.15.
%
% Inputs:
%   x_fdoa      (Ndim x N) array of FDOA sensor positions
%   v_fdoa      (Ndim x N) array of FDOA sensor velocities
%   xs          (Ndim x M) array of source positions over which to 
%               calculate CRLB
%   C           Covariance matrix for range rate estimates at the N
%               FDOA sensors [(m/s)^2]
%   a_grad      Function handle that returns the gradient of the
%               equality constraints a(x) as an nDim x nConstraint 
%               matrix
%   ref_idx     Scalar index of reference sensor, or nDim x nPair
%               matrix of sensor pairings
%   
% Outputs:
%   crlb    Lower bound on the error covariance matrix for an unbiased
%           FDOA estimator (Ndim x Ndim)
%
% Nicholas O'Donoughue
% 17 November 2021

% Parse inputs
if nargin < 5 || ~exist('ref_idx','var')
    ref_idx = [];
end
[n_dim, n_sensor] = size(x_fdoa);
n_source = size(xs,2);

% Define the Jacobian function
J = @(x) fdoa.jacobian(x_fdoa,v_fdoa,x,ref_idx);

% Resample covariance matrix
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

% Initialize output variable
crlb = zeros([n_dim,n_dim,n_source]);

warning('off','MATLAB:nearlySingularMatrix');

% Repeat CRLB for each of the n_source test positions
for idx =1:n_source
    this_x = xs(:,idx);
    
    % Evaluate the Jacobian
    J_i = J(this_x);
    
    % Evaluate the Gradient of the Constraint matrix
    A_i = a_grad(this_x);

    % Compute the Fisher Information Matrix
    if do_decomp
        F = J_i/C_d*J_i'; % Ndim x Ndim
    else
        F = J_i*C_inv*J_i';
    end
    
    if any(isnan(F(:))) || any(isinf(F(:)))
        % Problem is ill defined, Fisher Information Matrix cannot be
        % inverted
        crlb(:,:,idx) = NaN;
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