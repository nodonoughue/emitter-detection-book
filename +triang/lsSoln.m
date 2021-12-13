function [x,x_full] = lsSoln(x_aoa,psi,C,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress)
% [x,x_full] = lsSoln(x_aoa,psi,C,x_init,epsilon,max_num_iterations,...
%                                           force_full_calc,plot_progress)
%
% Computes the least square solution given a set of angle of arrival
% estimates.
%
% Inputs:
%   x_aoa               Sensor positions [m]
%   psi                 Measurement vector [radians]
%   C                   AOA measurement error covariance [radians^2]
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
% 1 July 2019

% Parse inputs
if nargin < 8
    plot_progress = [];
end

if nargin < 7
    force_full_calc = [];
end

if nargin < 6
    max_num_iterations = [];
end

if nargin < 5
    epsilon = [];
end

% Set up function handles
y = @(x) psi - triang.measurement(x_aoa, x);
J = @(x) triang.jacobian(x_aoa, x);

% Call the generic Least Square solver
[x,x_full] = utils.lsSoln(y,J,C,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress);