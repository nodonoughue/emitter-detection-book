function [x,x_full] = lsSoln(x_aoa,x_tdoa,x_fdoa,v_fdoa,z,C,x_init,epsilon,max_num_iterations, force_full_calc, plot_progress)
% [x,x_full] = lsSoln(x0,rho,C,xs_init,epsilon,numIterations,
%                                      force_full_calc, plot_progress)
%
% Computes the least square solution for combined AOA, TDOA, and
% FDOA processing.
%
% Inputs:
%   
%   x0                  Sensor positions [m]
%   rho                 Combined measurement vector
%   C                   Combined Error Covariance Matrix
%   xs_init             Initial estimate of source position [m]
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

% Initialize measurement error and Jacobian function handles
y = @(x) z - hybrid.measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa, x);
J = @(x) hybrid.jacobian(x_aoa, x_tdoa, x_fdoa, v_fdoa, x);

% Call the generic Least Square solver
[x,x_full] = utils.lsSoln(y,J,C,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress);