function [x_est,A,x_grid] = mlSoln(x_tx,x_rx,rho,C,x_ctr,search_size,epsilon,ref_idx)
% function [x_est,A,x_grid] = mlSoln(x_tx,x_rx,rho,C,x_ctr,search_size,...
%                                                       epsilon,ref_idx)
%
% Computes the Maximum Likelihood estimate for source position from a set
% of PCL bistatic range measurements.
%
% TODO: Incorporate range rate measurements
% TODO: Incorporate angle of arrival measurements
%
% Inputs:
%   
%   x_tx                Transmitter positions [m]
%   x_rx                Receiver positions [m]
%   rho                 Bistatic Range Measurements [m]
%   C                   Measurement Error Covariance Matrix [m^2]
%   x_ctr               Center of search grid [m]
%   search_size         2-D vector of search grid sizes [m]
%   epsilon             Desired resolution of search grid [m]
%   ref_idx             Matrix of tx/rx pairing indices (in the order
%                       they're used in C).  If ignored, then all pairwise
%                       measurements are used (nTx x nRx)
%
% OUTPUTS:
%   x_est           Estimated source position [m]
%   A               Likelihood computed across the entire set of candidate
%                   source positions
%   x_grid          Candidate source positions
%
% Nicholas O'Donoughue
% 10 September 2021

if nargin < 8 || ~exist('ref_idx','var')
    ref_idx = [];
end

% Set up function handle
ell = @(x) pcl.loglikelihood(x_tx,x_rx,rho,C,x,ref_idx);

% Call the util function
[x_est,A,x_grid] = utils.mlSoln(ell,x_ctr,search_size,epsilon);