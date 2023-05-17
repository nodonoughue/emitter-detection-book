function crlb = computeCRLBpriors(x_aoa,xs,C,fim_prior,do2DAoA)
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
if nargin < 5 || isempty(do2DAoA)
    do2DAoA = ~(size(C,1)==size(x_aoa,2)); % If the cov mtx is 1 per sensor, then it's not a 2D AOA problem 
end

% Set up Jacobian function
J = @(x) triang.jacobian(x_aoa,x,do2DAoA);

% Call the generic solver with the prior information
crlb = utils.computeCRLB(xs, C_tilde, J, fim_prior);