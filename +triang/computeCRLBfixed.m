function crlb = computeCRLBfixed(x_aoa,xs,C,a_grad,do2DAoA)
% crlb = computeCRLBfixed(x_aoa,xs,C,a_grad,do2DAoA)
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
%   do2DAoA         Boolean flag indicating whether 2D AOA (azimuth and
%                   elevation) is to be used
%
% Outputs:
%   crlb    Lower bound on the error covariance matrix for an unbiased
%           AOA estimator (Ndim x Ndim)
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 5 || isempty(do2DAoA)
    do2DAoA = ~(size(C,1)==size(x_aoa,2)); % If the cov mtx is 1 per sensor, then it's not a 2D AOA problem 
end

% Set up Jacobian function
J = @(x) triang.jacobian(x_aoa,x,do2DAoA);

% Call the generic solver with the equality constraints
crlb = utils.computeCRLB(xs, C, J, [], a_grad);