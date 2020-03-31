function [x_est,A,x_grid] = mlSoln(x_sensor,psi,C,x_ctr,search_size,epsilon)
% [x_est,A,x_grid] = mlSoln(x_sensor,psi,C,x_ctr,search_size,epsilon)
%
% Construct the ML Estimate by systematically evaluating the log
% likelihood function at a series of coordinates, and returning the index
% of the maximum.  Optionally returns the full set of evaluated
% coordinates, as well.
%
% INPUTS:
%   x_sensor    Sensor positions [m]
%   psi         Measurement vector [radians]
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

% Set up function handle
ell = @(x) triang.loglikelihood(x_sensor, psi,C,x);

% Call the util function
[x_est,A,x_grid] = utils.mlSoln(ell,x_ctr,search_size,epsilon);