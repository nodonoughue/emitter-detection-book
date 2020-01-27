function [x,x_full] = gdSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, z,C,x_init,alpha,beta,epsilon,max_num_iterations,force_full_calc,plot_progress)
% [x,x_full] = gdSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, z,C,x_init,alpha,...
%            beta,epsilon,max_num_iterations,force_full_calc,plot_progress)
%
% Computes the gradient descent solution for TDOA processing.
%
% Inputs:   
%   x_aoa               AOA sensor positions
%   x_tdoa              TDOA sensor positions
%   x_fdoa              FDOA sensor positions
%   v_fdoa              FDOA sensor velocities
%   z                   Measurement vector
%   C                   Combined error covariance matrix
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

y = @(x) z - hybrid.measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa, x);
J = @(x) hybrid.jacobian(x_aoa, x_tdoa, x_fdoa, v_fdoa, x);

[x,x_full] = utils.gdSoln(y,J,C,x_init,alpha,beta,epsilon,max_num_iterations,force_full_calc,plot_progress);
