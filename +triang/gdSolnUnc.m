function [x,x_full,alpha_est,beta_est] = gdSolnUnc(x_aoa, z,C,x_init,a,b,epsilon,max_num_iterations,force_full_calc,plot_progress)
% [x,x_full,alpha_est,beta_est] = gdSolnUnc(x_aoa, z,C,x_init,a,
%                        b,epsilon,max_num_iterations,force_full_calc,
%                        plot_progress)
%
% Computes the gradient descent solution for AOA processing.
%
% Inputs:   
%   x_aoa               AOA sensor positions
%   z                   Measurement vector
%   C                   Combined error covariance matrix
%   x_init              Initial estimate of source position [m]
%   a                   Backtracking line search parameter
%   b                   Backtracking line search parameter
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
%   alpha_est       Array with bias estimates
%   beta_est        Array with estimated sensor positions
%
% Nicholas O'Donoughue
% 17 Feb 2022

% Parse inputs
n_dim = size(x_aoa,1);
n_aoa = size(x_aoa,2);

assert(size(C,1) == n_aoa || size(C,1) == 2*n_aoa,'Unable to determine if AOA measurements are 1D or 2D');
do2DAoA = size(C,1) == 2*n_aoa;
if do2DAoA
    m_aoa = 2*n_aoa;
else
    m_aoa = n_aoa;
end

% Initialize measurement error and Jacobian function handles
% theta vector contains x, alpha, and beta.  Let's define the
% indices
x_ind = 1:n_dim;
alpha_ind = x_ind(end) + (1:m_aoa);
beta_ind = alpha_ind(end) + (1:n_dim*n_aoa);

y = @(theta) z - triang.measurement(reshape(theta(beta_ind),n_dim,n_aoa), ...   % x_aoa
                                    reshape(theta(x_ind),n_dim,1), ...            % x
                                    do2DAoA, ...
                                    reshape(theta(alpha_ind),m_aoa,1));      % alpha_aoa

J = @(theta) triang.jacobianUnc(reshape(theta(beta_ind),n_dim,n_aoa), ...   % x_aoa
                                reshape(theta(x_ind),n_dim,1), ...            % x
                                do2DAoA, ...
                                reshape(theta(alpha_ind),m_aoa,1));      % alpha_aoa
                                
% Build the initial theta vector
th_init = [x_init; zeros(m_aoa,1); x_aoa(:)];

% Call the generic Least Square solver
[th,th_full] = utils.gdSoln(y,J,C,th_init,a,b,epsilon,max_num_iterations,force_full_calc,plot_progress);

% Grab the x coordinates
x = th(x_ind);
x_full = th_full(x_ind,:);
alpha_est = th(alpha_ind);
beta_est = th(beta_ind);
