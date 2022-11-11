function crlb = computeCRLBpriors(x_aoa,x0,C,fim_prior,do2DAoA)
% crlb = computeCRLBpriors(x_aoa,xs,C,fim_prior,do2DAoA)
%
% Computes the CRLB on position accuracy for source at location xs and
% sensors at locations in x_aoa (Ndim x N).  C is an NxN matrix of TOA
% covariances at each of the N sensors.
%
% This version computes the CRLB in the presence of a statistical prior,
% created using the util.makePrior function, which returns not only a
% function handle for the prior, but the FIM for the log of the prior. The
% combined CRLB is computed according to (5.38).
%
% Inputs:
%   x_aoa           (Ndim x N) array of AOA sensor positions
%   x0              (Ndim x M) array of source positions over which to 
%                   calculate CRLB
%   C               AOA covariance matrix (radians^2)
%   fim_prior   Function handle that computes the Fisher Information Matrix
%               for the statistical prior
%   do2DAoA         Boolean flag indicating whether 2D AOA (az/el) should
%                   be assumed
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

if nargin < 5 || isempty(do2DAoA)
    do2DAoA = ~(size(C,1)==size(x_aoa,2)); % If the cov mtx is 1 per sensor, then it's not a 2D AOA problem 
end

% Set up Jacobian function
J = @(x) triang.jacobian(x_aoa,x,do2DAoA);

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
    
    % Compute the Fisher Information Matrix of the prior and add to F
    F = F + fim_prior(this_x);
    
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
