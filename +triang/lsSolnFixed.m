function [x,x_full] = lsSolnFixed(x_aoa,psi,C,x_init,a,tol,epsilon,max_num_iterations,force_full_calc,plot_progress)
% [x,x_full] = lsSolnFixed(x_aoa,psi,C,x_init,a,tol,epsilon,max_num_iterations,...
%                                           force_full_calc,plot_progress)
%
% Computes the least square solution given a set of angle of arrival
% estimates.
%
% Utilized the utils.constraints package to accept equality constraints (a)
% with tolerance (tol).
%
% Inputs:
%   x_aoa               Sensor positions [m]
%   psi                 Measurement vector [radians]
%   C                   AOA measurement error covariance [radians^2]
%   x_init              Initial source position estimate [m]
%   a                   Array of equality constraint function handles
%   tol                 Tolerance for equality constraints
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
% 14 November 2021

% Parse inputs
if nargin < 10
    plot_progress = [];
end

if nargin < 9
    force_full_calc = [];
end

if nargin < 8
    max_num_iterations = [];
end

if nargin < 7
    epsilon = [];
end

if nargin < 6
    tol = [];
end

% Set up function handles
y = @(x) psi - triang.measurement(x_aoa, x);
J = @(x) triang.jacobian(x_aoa, x);

% Call the generic Least Square solver
[x,x_full] = utils.constraints.lsSolnFixed(y,J,C,x_init,a,tol,...
    epsilon,max_num_iterations,force_full_calc,plot_progress);