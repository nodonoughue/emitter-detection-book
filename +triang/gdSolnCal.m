function [x,x_full,alpha_est,beta_est] = gdSolnCal(x_aoa,z,x_cal, z_cal, C,x_init,a,b,epsilon,max_num_iterations,force_full_calc,plot_progress)
% [x,x_full,alpha_est,beta_est] = gdSolnCal(x_aoa,z,x_cal,z_cal,C,x_init,
%               a,b,epsilon,max_num_iterations,force_full_calc,plot_progress)
%
% Computes the gradient descent solution for AOA processing using 
% calibration measurements to first estimate the unknown sensor 
% measurement biases (alpha_est) and sensor positions (beta_est).
%
% Inputs:   
%   x_aoa               AOA sensor positions
%   z                   Measurement vector
%   x_cal               Calibration emitter positions (n_dim x n_cal)
%   z_cal               Calibration emitter measurements (n_msmt x n_cal)
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
%
% Nicholas O'Donoughue
% 17 Feb 2022

% Parse inputs

% Parse inputs sizes
n_dim = size(x_aoa,1);
n_aoa = size(x_aoa,2);
n_cal = size(x_cal,2);

assert(size(C,1) == n_aoa || size(C,1) == 2*n_aoa ,'Unable to determine if AOA measurements are 1D or 2D');
do2DAoA = size(C,1) == 2*n_aoa;
if do2DAoA
    m_aoa = 2*n_aoa;
else
    m_aoa = n_aoa;
end

% Expand the noise for cal measurements (n_cal independent measurements)
C_cal = kron(eye(n_cal),C);

%% Estimate Measurement Bias
%  Assume provided sensor positions are accurate
%  Use known calibration emitter positions
y_a = @(alpha) reshape(z_cal - triang.measurement(x_aoa,x_cal,do2DAoA,alpha),...
                                                  m_aoa*n_cal,1);
    % Response will be n_msmt x n_cal, reshape to n_msmt*n_cal x 1

J_a = @(alpha) reshape(triang.grad_a(x_aoa, x_cal, do2DAoA, alpha), ...
                                     m_aoa,m_aoa*n_cal);
    % Response will be n_alpha x n_msmt x n_cal, reshape to n_dim x n_msmt*n_cal

% Build the initial theta vector
alpha_init = zeros(m_aoa,1);

% Call the generic Least Square solver
alpha_est = utils.gdSoln(y_a,J_a,C_cal,alpha_init,a,b,epsilon,max_num_iterations,force_full_calc,plot_progress);

%% Estimate Sensor Position Error
%  Use estimated alpha, and known cal positions
%  Estimate sensor positions
y_b = @(beta) reshape(z_cal - triang.measurement(reshape(beta,n_dim,n_aoa), ...   % x_aoa
                                                 x_cal, do2DAoA, alpha_est), ...
                                                 m_aoa*n_cal,1);

J_b = @(beta) reshape(triang.grad_b(reshape(beta,n_dim,n_aoa), ...   % x_aoa
                                    x_cal, do2DAoA, alpha_est),...
                                    n_dim*n_aoa,m_aoa*n_cal);

% Build the initial vector
beta_init = x_aoa(:);

beta_est = utils.gdSoln(y_b,J_b,C_cal,beta_init,a,b,epsilon,max_num_iterations,force_full_calc,plot_progress);
x_aoa_est = reshape(beta_est,n_dim,n_aoa);

%% Estimate Target Position
%  Use alpha_est and beta_est with tgt measurements z to estimate

y = @(x) z - triang.measurement(x_aoa_est, x, do2DAoA, alpha_est);
J = @(x) triang.grad_x(x_aoa_est, x, do2DAoA, alpha_est);

[x,x_full] = utils.gdSoln(y,J,C_tilde,x_init,a,b,epsilon,max_num_iterations,force_full_calc,plot_progress);