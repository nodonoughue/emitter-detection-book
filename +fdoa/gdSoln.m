function [x,x_full] = gdSoln(x_fdoa,v_fdoa,rho_dot,C,x_init,alpha,beta,epsilon,max_num_iterations,force_full_calc,plot_progress)
% [x,x_full] = gradDescentSoln(x_fdoa,v_fdoa,rho_dot,C,x_init,alpha,...
%           beta,epsilon,max_num_iterations,force_full_calc,plot_progress)
%
% Computes the gradient descent solution for FDOA processing.
%
% Inputs:
%   
%   x_fdoa              FDOA sensor positions [m]
%   v_fdoa              FDOA sensor velocities [m/s]
%   rho_dot             Measurement vector
%   C                   FDOA error covariance matrix
%   x_init              Initial estimate of source position [m]
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
% 1 July 2019

% Initialize measurement error and jacobian functions
y = @(x) rho_dot- fdoa.measurement(x_fdoa, v_fdoa, x);
J = @(x) fdoa.jacobian(x_fdoa,v_fdoa,x);

% Call generic Gradient Descent solver
[x,x_full] = utils.gdSoln(y,J,C,x_init,alpha,beta,epsilon,max_num_iterations,force_full_calc,plot_progress);
