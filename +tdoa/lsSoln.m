function [x,x_full] = lsSoln(x_tdoa,rho,C,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress,ref_idx)
% [x,x_full] = lsSoln(x_tdoa,rho,C,x_init,epsilon,max_num_iterations,...
%                               force_full_calc, plot_progress,ref_idx)
%
% Computes the least square solution for TDOA processing.
%
% Inputs:
%   
%   x_tdoa              Sensor positions [m]
%   rho                 Range-Difference Measurements [m]
%   C                   Range-Difference Error Covariance Matrix [m^2]
%   x_init              Initial source position estimate [m]
%   epsilon             Desired estimate resolution [m]
%   max_num_iterations  Maximum number of iterations to perform
%   force_full_calc     Boolean flag to force all iterations (up to
%                       max_num_iterations) to be computed, regardless
%                       of convergence (DEFAULT = False)
%   plot_progress       Boolean flag dictacting whether to plot
%                       intermediate solutions as they are derived 
%                       (DEFAULT = False).
%   ref_idx             Scalar index of reference sensor, or nDim x nPair
%                       matrix of sensor pairings           
% 
% Outputs:
%   x               Estimated source position
%   x_full          Iteration-by-iteration estimated source positions
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 9 || ~exist('ref_idx','var')
    ref_idx = [];
end

if nargin < 8 || ~exist('plot_progress','var')
    plot_progress = false;
end

if nargin < 7 || ~exist('force_full_calc','var')
    force_full_calc = false;
end

if nargin < 6 || ~exist('max_num_iterations','var')
    max_num_iterations = [];
end

if nargin < 5 || ~exist('epsilon','var')
    epsilon = [];
end

% Initialize measurement error and Jacobian function handles
y = @(x) rho- tdoa.measurement(x_tdoa, x, ref_idx);
J = @(x) tdoa.jacobian(x_tdoa, x, ref_idx);

% Call the generic Least Square solver
[x,x_full] = utils.lsSoln(y,J,C,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress);