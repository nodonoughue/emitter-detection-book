function [x_est,A,x_grid] = mlSolnPrior(x_fdoa,v_fdoa,zeta,C,prior,x_ctr,search_size,epsilon,lambda,fdoa_ref_idx)
% [x_est,A,x_grid] = mlSolnPrior(x_fdoa,v_fdoa,zeta,C,prior,x_ctr,...
%                               search_size,epsilon,lambda,fdoa_ref_idx)
%
% Construct the ML Estimate by systematically evaluating the log
% likelihood function at a series of coordinates, and returning the index
% of the maximum.  Optionally returns the full set of evaluated
% coordinates, as well.
%
% Accepts a prior distribution on the position of the target x.
%
% INPUTS:
%   x_fdoa          FDOA sensor positions [m]
%   v_fdoa          FDOA sensor velocities [m/s]
%   zeta            Combined measurement vector
%   C               Combined measurement error covariance matrix
%   prior           Function handle to prior distribution
%   x_ctr           Center of search grid [m]
%   search_size     2-D vector of search grid sizes [m]
%   epsilon         Desired resolution of search grid [m]
%   lambda          Optional weight factor for prior (0-1, default=.5)
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA
%
% OUTPUTS:
%   x_est           Estimated source position [m]
%   A               Likelihood computed across the entire set of candidate
%                   source positions
%   x_grid          Candidate source positions
%
% Nicholas O'Donoughue
% 24 November 2021

% Parse inputs
if nargin < 9 || ~exist('lambda','var')
    lambda = .5;
end
assert(isscalar(lambda) && isfinite(lambda) && lambda >= 0 && lambda <= 1,...
       'The prior distribution weight lambda must be a scalar between 0 and 1.');

if nargin < 10 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

% Set up function handle
ell = @(x) fdoa.loglikelihood(x_fdoa,v_fdoa,zeta,C,x,fdoa_ref_idx);

% Adjust likelihood for prior
ell_prior = @(x) (1-lambda)*ell(x) + lambda*log10(prior(x));

% Call the util function
[x_est,A,x_grid] = utils.mlSoln(ell_prior,x_ctr,search_size,epsilon);
