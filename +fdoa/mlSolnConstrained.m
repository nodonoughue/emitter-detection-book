function [x_est,A,x_grid] = mlSolnConstrained(x_fdoa,v_fdoa,rho_dot,C,x_ctr,search_size,epsilon,a,b,tol, ref_idx)
% [x_est,A,x_grid] = mlSolnConstrained(x_fdoa,v_fdoa,rho_dot,C,x_ctr,search_size,epsilon,a,b,tol, ref_idx)
%
% Construct the ML Estimate by systematically evaluating the log
% likelihood function at a series of coordinates, and returning the index
% of the maximum.  Optionally returns the full set of evaluated
% coordinates, as well.
%
% This version of mlSoln accepts a constraint function a(x) that must equal
% zero for all acceptable points, and a vector of function handles b that
% must all be <=0 zero. The cost returned for all points that
% violate these constraints is Inf.
%
% INPUTS:
%   x_fdoa      Sensor positions [m]
%   v_fdoa      Sensor velocities [m/s]
%   rho_dot     Measurement vector [Hz]
%   C           Measurement error covariance matrix
%   x_ctr       Center of search grid [m]
%   search_size 2-D vector of search grid sizes [m]
%   epsilon     Desired resolution of search grid [m]
%   ref_idx     Scalar index of reference sensor, or nDim x nPair
%               matrix of sensor pairings
%   a           Equality constraint function handle
%   b           Array of inequality constraint function handles
%   tol         Tolerance for equality constraint
%
% OUTPUTS:
%   x_est           Estimated source position [m]
%   A               Likelihood computed across the entire set of candidate
%                   source positions
%   x_grid          Candidate source positions
%
% Nicholas O'Donoughue
% 5 September 2021

if nargin < 11 || ~exist('ref_idx','var')
    ref_idx = [];
end

if nargin < 10 || ~exist('tol','var')
    tol = [];
end

if nargin < 9 || ~exist('b','var')
    b = [];
end

if nargin < 8 || ~exist('a','var')
    a = [];
end

% Set up function handle
ell = @(x) fdoa.loglikelihood(x_fdoa,v_fdoa,rho_dot,C,x,ref_idx);

% Call the util function
[x_est,A,x_grid] = utils.constraints.mlSolnConstrained(ell,x_ctr,search_size,epsilon,a,b,tol);