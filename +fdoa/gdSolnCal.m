function [x,x_full,alpha_est,beta_est] = gdSolnCal(x_fdoa, v_fdoa, z,x_cal, z_cal, C,x_init,a,b,epsilon,max_num_iterations,force_full_calc,plot_progress,fdoa_ref_idx)
% [x,x_full,alpha_est,beta_est] = gdSolnCal(x_fdoa, v_fdoa, 
%                   z,x_cal,z_cal,C,x_init,a,b,epsilon,max_num_iterations,
%                   force_full_calc,plot_progress,fdoa_ref_idx)
%
% Computes the gradient descent solution for FDOA processing using 
% calibration measurements to first estimate the unknown sensor 
% measurement biases (alpha_est) and sensor positions (beta_est).
%
% Inputs:   
%   x_fdoa              FDOA sensor positions
%   v_fdoa              FDOA sensor velocities
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
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA
%
% Outputs:
%   x               Estimated source position
%   x_full          Iteration-by-iteration estimated source positions
%
% Nicholas O'Donoughue
% 17 Feb 2022

% Parse inputs
if nargin < 14 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

% Parse inputs sizes
n_dim = size(x_fdoa,1);
n_fdoa = size(x_fdoa,2);
n_cal = size(x_cal,2);

[test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(fdoa_ref_idx,n_fdoa);
m_fdoa = numel(test_idx_vec);

%% Re-sample Noise

% For now, we assume the AOA is independent of TDOA/FDOA
C_tilde = utils.resampleCovMtx(C, test_idx_vec, ref_idx_vec);

% Expand the noise for cal measurements (n_cal independent measurements)
C_cal = kron(eye(n_cal),C_tilde);

%% Estimate Measurement Bias
%  Assume provided sensor positions are accurate
%  Use known calibration emitter positions
y_a = @(alpha) reshape(z_cal - fdoa.measurement(x_fdoa,v_fdoa,...
                                                x_cal,fdoa_ref_idx,alpha),...
                                                m_fdoa*n_cal,1);
    % Response will be n_msmt x n_cal, reshape to n_msmt*n_cal x 1

J_a = @(alpha) reshape(fdoa.grad_a(x_fdoa, v_fdoa, x_cal, ...
                                     fdoa_ref_idx, alpha),...
                                     n_fdoa,m_fdoa*n_cal);
    % Response will be n_alpha x n_msmt x n_cal, reshape to n_dim x n_msmt*n_cal

% Build the initial theta vector
alpha_init = zeros(n_fdoa,1);

% Call the generic Least Square solver
alpha_est = utils.gdSoln(y_a,J_a,C_cal,alpha_init,a,b,epsilon,max_num_iterations,force_full_calc,plot_progress);

%% Estimate Sensor Position Error
%  Use estimated alpha, and known cal positions
%  Estimate sensor positions
beta_fx_ind = 1:n_dim*n_fdoa;
beta_fv_ind = beta_fx_ind(end)+ (1:n_dim*n_fdoa);
y_b = @(beta) reshape(z_cal - fdoa.measurement(reshape(beta(beta_fx_ind),n_dim,n_fdoa), ... % x_fdoa
                                               reshape(beta(beta_fv_ind),n_dim,n_fdoa), ... % v_fdoa
                                               x_cal, ...            % x
                                               fdoa_ref_idx, alpha_est), ...
                                               m_fdoa*n_cal,1);

J_b = @(beta) reshape(fdoa.grad_b(reshape(beta(beta_fx_ind),n_dim,n_fdoa), ... % x_fdoa
                                  reshape(beta(beta_fv_ind),n_dim,n_fdoa), ... % v_fdoa
                                  x_cal, ...            % x
                                  fdoa_ref_idx, alpha_est), ...
                                  n_dim*2*n_fdoa,m_fdoa*n_cal);

% Build the initial vector
beta_init = [x_fdoa(:); v_fdoa(:)];

beta_est = utils.gdSoln(y_b,J_b,C_cal,beta_init,a,b,epsilon,max_num_iterations,force_full_calc,plot_progress);
x_fdoa_est = reshape(beta_est(beta_fx_ind),n_dim,n_fdoa);
v_fdoa_est = reshape(beta_est(beta_fv_ind),n_dim,n_fdoa);

%% Estimate Target Position
%  Use alpha_est and beta_est with tgt measurements z to estimate

y = @(x) z - fdoa.measurement(x_fdoa_est, v_fdoa_est, x,fdoa_ref_idx, alpha_est);
J = @(x) fdoa.grad_x(x_fdoa_est, v_fdoa_est, x, fdoa_ref_idx, alpha_est);

[x,x_full] = utils.gdSoln(y,J,C_tilde,x_init,a,b,epsilon,max_num_iterations,force_full_calc,plot_progress);