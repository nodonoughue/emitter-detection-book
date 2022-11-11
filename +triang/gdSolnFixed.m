function [x,x_full] = gdSolnFixed(x_aoa,psi,C,x_init,a,tol,alpha,beta,epsilon,max_num_iterations,force_full_calc,plot_progress)
% [x,x_full] = gdSolnFixed(x_aoa,psi,C,x_init,a,tol,alpha,beta,epsilon,...
%                       max_num_iterations,force_full_calc,plot_progress)
%
% Computes the gradient descent solution for AOA processing.
%
% Utilized the utils.constraints package to accept equality constraints (a)
% with tolerance (tol).
%
% Inputs:
%   
%   x_aoa               AOA sensor positions [m]
%   psi                 Measurement vector
%   C                   Combined error covariance matrix
%   x_init              Initial estimate of source position [m]
%   a                   Array of equality constraint function handles
%   tol                 Tolerance for equality constraints
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
if nargin < 12
    plot_progress = [];
end

if nargin < 11
    force_full_calc = [];
end

if nargin < 10
    max_num_iterations = [];
end

if nargin < 9
    epsilon = [];
end

if nargin < 8
    beta = [];
end

if nargin < 7
    alpha = [];
end

if nargin < 6
    tol = [];
end

% Determine if measurement should be 2D (az/el) or 1D (az only)
if numel(psi) == size(x_aoa,2)
    do2DAoA = false;
else
    do2DAoA = true;
end

% Initialize measurement error and Jacobian functions
y = @(x) psi - triang.measurement(x_aoa, x,do2DAoA);
J = @(x) triang.jacobian(x_aoa,x,do2DAoA);

% Call the generic Gradient Descent solver
[x,x_full] = utils.constraints.gdSolnFixed(y,J,C,x_init,a,tol,...
    alpha,beta,epsilon,max_num_iterations,force_full_calc,plot_progress);
