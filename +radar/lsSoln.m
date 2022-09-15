function [x,x_full] = lsSoln(x_rdr,rho,C,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress)
% [x,x_full] = lsSoln(x_rdr,rho,C,x_init,epsilon,max_num_iterations,...
%                               force_full_calc, plot_progress,ref_idx)
%
% Computes the least square solution for radar processing -- range
% measurements only for now.
%
% TODO: Incorporate range rate measurements
% TODO: Incorporate angle of arrival measurements
%
% Inputs:
%   
%   x_rdr               Radar positions [m]
%   rho                 Measurements [m]
%   C                   Measurement Error Covariance Matrix [m^2]
%   x_init              Initial source position estimate [m]
%   epsilon             Desired estimate resolution [m]
%   max_num_iterations  Maximum number of iterations to perform
%   force_full_calc     Boolean flag to force all iterations (up to
%                       max_num_iterations) to be computed, regardless
%                       of convergence (DEFAULT = False)
%   plot_progress       Boolean flag dictacting whether to plot
%                       intermediate solutions as they are derived 
%                       (DEFAULT = False).
% 
% Outputs:
%   x               Estimated source position
%   x_full          Iteration-by-iteration estimated source positions
%
% Nicholas O'Donoughue
% 10 September 2021

% Parse inputs
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
y = @(x) rho - radar.measurement(x_rdr, x);
J = @(x) radar.jacobian(x_rdr, x);

% Call the generic Least Square solver
[x,x_full] = utils.lsSoln(y,J,C,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress);