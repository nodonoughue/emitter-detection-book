function [x,x_full] = gdSolnBounded(x_aoa,psi,C,x_init,b,alpha,beta,epsilon,max_num_iterations,force_full_calc,plot_progress)
% [x,x_full] = gdSolnBounded(x_aoa,psi,C,x_init,b,alpha,beta,epsilon,...
%                       max_num_iterations,force_full_calc,plot_progress)
%
% Computes the gradient descent solution for AOA processing.
%
% Utilized the utils.constraints package to accept inequality constraints 
% (b).
%
% Inputs:
%   
%   x_aoa               AOA sensor positions [m]
%   psi                 Measurement vector
%   C                   Combined error covariance matrix
%   x_init              Initial estimate of source position [m]
%   b                   Array of inequality constraint function handles
%   alpha               Backtracking line search parameter
%   beta                Backtracking line search parameter
%   epsilon             Desired position error tolerance (stopping 
%                       condition)
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
if nargin < 11
    plot_progress = [];
end

if nargin < 10
    force_full_calc = [];
end

if nargin < 9
    max_num_iterations = [];
end

if nargin < 8
    epsilon = [];
end

if nargin < 7
    beta = [];
end

if nargin < 6
    alpha = [];
end

% Initialize measurement error and Jacobian functions
y = @(x) psi - triang.measurement(x_aoa, x);
J = @(x) triang.jacobian(x_aoa,x);

% Call the generic Gradient Descent solver
[x,x_full] = utils.constraints.gdSolnBounded(y,J,C,x_init,b,...
    alpha,beta,epsilon,max_num_iterations,force_full_calc,plot_progress);
