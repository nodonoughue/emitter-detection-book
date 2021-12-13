function crlb = computeCRLB(x_aoa,x0,C)
% crlb = computeCRLB(x_aoa,xs,C)
%
% Computes the CRLB on position accuracy for source at location xs and
% sensors at locations in x_aoa (Ndim x N).  C is an NxN matrix of TOA
% covariances at each of the N sensors.
%
% Inputs:
%   x_aoa           (Ndim x N) array of AOA sensor positions
%   x0              (Ndim x M) array of source positions over which to 
%                   calculate CRLB
%   C               AOA covariance matrix (radians^2)
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

% Repeat CRLB for each of the n_source test positions
for idx =1:n_source
    this_x = x0(:,idx);
    
    % Compute Jacobian matrix
    J_i = J(this_x);
    
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
        crlb(:,:,idx) = pinv(F);
    end
end
