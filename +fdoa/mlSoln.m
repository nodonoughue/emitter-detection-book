function [x_est,A,x_grid] = mlSoln(x_fdoa,v_fdoa,rho_dot,C,x_ctr,search_size,epsilon,ref_idx)
% [x_est,A,x_grid] = mlSoln(x_fdoa,v_fdoa,rho_dot,C,x_ctr,search_size,...
%                                                       epsilon,ref_idx)
%
% Construct the ML Estimate by systematically evaluating the log
% likelihood function at a series of coordinates, and returning the index
% of the maximum.  Optionally returns the full set of evaluated
% coordinates, as well.
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
%
% OUTPUTS:
%   x_est           Estimated source position [m]
%   A               Likelihood computed across the entire set of candidate
%                   source positions
%   x_grid          Candidate source positions
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 8 || ~exist('ref_idx','var')
    ref_idx = [];
end

% Set up function handle
ell = @(x) fdoa.loglikelihood(x_fdoa,v_fdoa, rho_dot,C,x,ref_idx);

% Call the util function
[x_est,A,x_grid] = utils.mlSoln(ell,x_ctr,search_size,epsilon);