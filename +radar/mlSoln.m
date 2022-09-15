function [x_est,A,x_grid] = mlSoln(x_rdr,rho,C,x_ctr,search_size,epsilon)
% function [x_est,A,x_grid] = mlSoln(x_rdr,rho,C,x_ctr,search_size,epsilon)
%
% Computes the Maximum Likelihood estimate for source position from a set
% of monostatic range measurements.
%
% TODO: Incorporate range rate measurements
% TODO: Incorporate angle of arrival measurements
%
% Inputs:
%   
%   x_rdr               Radar positions [m]
%   rho                 Bistatic Range Measurements [m]
%   C                   Measurement Error Covariance Matrix [m^2]
%   x_ctr               Center of search grid [m]
%   search_size         2-D vector of search grid sizes [m]
%   epsilon             Desired resolution of search grid [m]
%
% OUTPUTS:
%   x_est           Estimated source position [m]
%   A               Likelihood computed across the entire set of candidate
%                   source positions
%   x_grid          Candidate source positions
%
% Nicholas O'Donoughue
% 15 September 2022

% Set up function handle
ell = @(x) radar.loglikelihood(x_rdr,[],rho,C,x,[],0);

% Call the util function
[x_est,A,x_grid] = utils.mlSoln(ell,x_ctr,search_size,epsilon);