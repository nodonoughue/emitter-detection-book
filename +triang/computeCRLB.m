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
%           TDOA estimator (Ndim x Ndim)
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
n_dim = size(x_aoa,1);
n_source = size(x0,2);

% Set up Jacobian function
J = @(x) triang.jacobian(x_aoa,x);

% Preprocess covariance matrix
C_d = decomposition(C);

% Initialize output variable
crlb = zeros([n_dim,n_dim,n_source]);

% Repeat CRLB for each of the n_source test positions
for idx =1:n_source
    this_x = x0(:,idx);
    J_i = J(this_x);
    F = J_i/C_d*J_i'; % Ndim x Ndim
    
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
