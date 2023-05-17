function crlb = computeCRLB(x_aoa,xs,C,do2DAoA)
% crlb = computeCRLB(x_aoa,xs,C,do2DAoA)
%
% Computes the CRLB on position accuracy for source at location xs and
% sensors at locations in x_aoa (Ndim x N).  C is an NxN matrix of TOA
% covariances at each of the N sensors.
%
% Inputs:
%   x_aoa           (Ndim x N) array of AOA sensor positions
%   xs              (Ndim x M) array of source positions over which to 
%                   calculate CRLB
%   C               AOA covariance matrix (radians^2)
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
if nargin < 4 || isempty(do2DAoA)
    do2DAoA = ~(size(C,1)==size(x_aoa,2)); % If the cov mtx is 1 per sensor, then it's not a 2D AOA problem 
end

% Set up Jacobian function
J = @(x) triang.jacobian(x_aoa,x,do2DAoA);

% Call Generic Solver
crlb = utils.computeCRLB(xs, C, J);