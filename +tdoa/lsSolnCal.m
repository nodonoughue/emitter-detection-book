function [x,x_full,alpha_est,beta_est] = lsSolnCal(x_tdoa, z,x_cal, z_cal, C,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress,tdoa_ref_idx)
% [x,x_full,alpha_est,beta_est] = lsSolnCal(x_tdoa, z,x_cal, z_cal, C,
%                           x_init,epsilon,max_num_iterations,
%                           force_full_calc,plot_progress,tdoa_ref_idx)
%
% Computes the least square solution for hybrid AOA, TDOA, and
% FDOA processing using calibration measurements to first estimate the
% unknown sensor measurement biases (alpha_est) and sensor positions
% (beta_est).
%
% Inputs:   
%   x_tdoa              TDOA sensor positions
%   z                   Measurement vector
%   x_cal               Calibration emitter positions (n_dim x n_cal)
%   z_cal               Calibration emitter measurements (n_msmt x n_cal)
%   C                   Combined error covariance matrix
%   x_init              Initial estimate of source position [m]
%   epsilon             Desired position error tolerance (stopping 
%                       condition)
%   max_num_iterations  Maximum number of iterations to perform
%   force_full_calc     Boolean flag to force all iterations (up to
%                       max_num_iterations) to be computed, regardless
%                       of convergence (DEFAULT = False)
%   plot_progress       Boolean flag dictacting whether to plot
%                       intermediate solutions as they are derived 
%                       (DEFAULT = False).
%   tdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for TDOA
%
% Outputs:
%   x               Estimated source position
%   x_full          Iteration-by-iteration estimated source positions
%
% Nicholas O'Donoughue
% 17 Feb 2022

% Parse inputs
if nargin < 11 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

% Parse inputs sizes
n_dim = size(x_tdoa,1);
n_tdoa = size(x_tdoa,2);
n_cal = size(x_cal,2);


[test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(tdoa_ref_idx,n_tdoa);
m_tdoa = numel(test_idx_vec);

%% Re-sample Noise
C_tilde = utils.resampleCovMtx(C, test_idx_vec, ref_idx_vec);

% Expand the noise for cal measurements (n_cal independent measurements)
C_cal = kron(eye(n_cal),C_tilde);

%% Estimate Measurement Bias
%  Assume provided sensor positions are accurate
%  Use known calibration emitter positions
y_a = @(alpha) reshape(z_cal - tdoa.measurement(x_tdoa,x_cal,tdoa_ref_idx,alpha),...
                                               (m_aoa+m_tdoa+m_fdoa)*n_cal,1);
    % Response will be n_msmt x n_cal, reshape to n_msmt*n_cal x 1

J_a = @(alpha) reshape(tdoa.grad_a(x_tdoa, x_cal, tdoa_ref_idx, alpha),...
                                   n_tdoa,m_tdoa*n_cal);
    % Response will be n_alpha x n_msmt x n_cal, reshape to n_dim x n_msmt*n_cal

% Build the initial theta vector
alpha_init = zeros(n_tdoa,1);

% Call the generic Least Square solver
alpha_est = utils.lsSoln(y_a,J_a,C_cal,alpha_init,epsilon,max_num_iterations,force_full_calc,plot_progress);

%% Estimate Sensor Position Error
%  Use estimated alpha, and known cal positions
%  Estimate sensor positions
y_b = @(beta) reshape(z_cal - tdoa.measurement(reshape(beta,n_dim,n_tdoa), ...  % x_tdoa
                                               x_cal, tdoa_ref_idx, ...
                                               alpha_est), ...
                                               m_tdoa*n_cal,1);

J_b = @(beta) reshape(tdoa.grad_b(reshape(beta,n_dim,n_tdoa), ...  % x_tdoa
                                  x_cal, tdoa_ref_idx, ...
                                  alpha_est), ...
                                  n_dim*n_tdoa,m_tdoa*n_cal);

% Build the initial vector
beta_init = x_tdoa(:);

beta_est = utils.lsSoln(y_b,J_b,C_cal,beta_init,epsilon,max_num_iterations,force_full_calc,plot_progress);
x_tdoa_est = reshape(beta_est,n_dim,n_tdoa);

%% Estimate Target Position
%  Use alpha_est and beta_est with tgt measurements z to estimate

y = @(x) z - tdoa.measurement(x_tdoa_est, x, tdoa_ref_idx, alpha_est);
J = @(x) tdoa.grad_x(x_tdoa_est, x, tdoa_ref_idx, alpha_est);

[x,x_full] = utils.lsSoln(y,J,C_tilde,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress);