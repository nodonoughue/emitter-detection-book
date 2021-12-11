function [x_est,A,x_grid] = mlSoln(x_tdoa,rho,C,x_ctr,search_size,epsilon,ref_idx)
% function [x_est,A,x_grid] = mlSoln(x_tdoa,rho,C,x_ctr,search_size,...
%                                                       epsilon,ref_idx)
%
% Computes the Maximum Likelihood estimate for source position from a set
% of TDOA measurements.
%
% INPUTS:
%   x_tdoa      Sensor positions [m]
%   rho         Measurement vectors [m]
%   C           Measurement error covariance matrix [m^2]
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

if nargin < 7 || ~exist('ref_idx','var')
    ref_idx = [];
end

% Set up function handle
ell = @(x) tdoa.loglikelihood(x_tdoa, rho,C,x,ref_idx);

% Call the util function
[x_est,A,x_grid] = utils.mlSoln(ell,x_ctr,search_size,epsilon);