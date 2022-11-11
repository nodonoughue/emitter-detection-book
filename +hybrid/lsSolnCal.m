function [x,x_full,alpha_est,beta_est] = lsSolnCal(x_aoa, x_tdoa, x_fdoa, v_fdoa, z,x_cal, z_cal, C,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress,tdoa_ref_idx,fdoa_ref_idx)
% [x,x_full,alpha_est,beta_est] = lsSolnCal(x_aoa, x_tdoa, x_fdoa, v_fdoa, 
%                   z,x_cal,z_cal,C,x_init,epsilon,max_num_iterations,
%                   force_full_calc,plot_progress,tdoa_ref_idx,fdoa_ref_idx)
%
% Computes the gradient descent solution for hybrid AOA, TDOA, and
% FDOA processing using calibration measurements to first estimate the
% unknown sensor measurement biases (alpha_est) and sensor positions
% (beta_est).
%
% Inputs:   
%   x_aoa               AOA sensor positions
%   x_tdoa              TDOA sensor positions
%   x_fdoa              FDOA sensor positions
%   v_fdoa              FDOA sensor velocities
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
if nargin < 14 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

if nargin < 15 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

% Parse inputs sizes
[n_dim_a,n_aoa] = size(x_aoa);
[n_dim_t,n_tdoa] = size(x_tdoa);
[n_dim_f,n_fdoa] = size(x_fdoa);
[n_dim_c,n_cal] = size(x_cal);

if n_aoa > 0
    n_dim = n_dim_a;
elseif n_tdoa > 0
    n_dim = n_dim_t;
elseif n_fdoa > 0
    n_dim = n_dim_f;
else
    error('At least one set of sensors required (AOA, TDOA, or FDOA)');
end
assert(n_dim_c == n_dim,'Error parsing input dimensions.');

assert(size(C,1) == n_aoa + n_tdoa + n_fdoa || size(C,1) == 2*n_aoa + n_tdoa + n_fdoa,'Unable to determine if AOA measurements are 1D or 2D');
do2DAoA = size(C,1) == 2*n_aoa;
if do2DAoA
    m_aoa = 2*n_aoa;
else
    m_aoa = n_aoa;
end

[tdoa_test_idx_vec, tdoa_ref_idx_vec] = utils.parseReferenceSensor(tdoa_ref_idx,n_tdoa);
[fdoa_test_idx_vec, fdoa_ref_idx_vec] = utils.parseReferenceSensor(fdoa_ref_idx,n_fdoa);

m_tdoa = numel(tdoa_test_idx_vec);
m_fdoa = numel(fdoa_test_idx_vec);

%% Re-sample Noise
% Parse the TDOA and FDOA reference indices together
test_idx_vec = cat(2,tdoa_test_idx_vec, n_tdoa + fdoa_test_idx_vec);
ref_idx_vec = cat(2,tdoa_ref_idx_vec, n_tdoa + fdoa_ref_idx_vec);

% For now, we assume the AOA is independent of TDOA/FDOA
C_aoa = C(1:n_aoa, 1:n_aoa);
C_tfdoa = C(n_aoa+1:end, n_aoa+1:end);
C_tilde = blkdiag(C_aoa, utils.resampleCovMtx(C_tfdoa, test_idx_vec, ref_idx_vec));

% Expand the noise for cal measurements (n_cal independent measurements)
C_cal = kron(eye(n_cal),C_tilde);

%% Estimate Measurement Bias
%  Assume provided sensor positions are accurate
%  Use known calibration emitter positions
alpha_a_ind = 1:m_aoa;
alpha_t_ind = m_aoa + 1:n_tdoa;
alpha_f_ind = m_aoa + (n_tdoa + 1:n_fdoa);

y_a = @(alpha) reshape(z_cal - hybrid.measurement(x_aoa,x_tdoa,x_fdoa,v_fdoa,...
                                                  x_cal,tdoa_ref_idx,fdoa_ref_idx,do2DAoA,...
                                                  alpha(alpha_a_ind),...
                                                  alpha(alpha_t_ind),...
                                                  alpha(alpha_f_ind)),...
                                                  (m_aoa+m_tdoa+m_fdoa)*n_cal,1);
    % Response will be n_msmt x n_cal, reshape to n_msmt*n_cal x 1

J_a = @(alpha) reshape(hybrid.grad_a(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_cal, ...
                                     tdoa_ref_idx, fdoa_ref_idx, do2DAoA, ...
                                     alpha(alpha_a_ind), ...
                                     alpha(alpha_t_ind), ...
                                     alpha(alpha_f_ind)),...
                                     m_aoa+n_tdoa+n_fdoa,(m_aoa+m_tdoa+m_fdoa)*n_cal);
    % Response will be n_alpha x n_msmt x n_cal, reshape to n_dim x n_msmt*n_cal

% Build the initial theta vector
alpha_init = zeros(m_aoa + n_tdoa + n_fdoa,1);

% Call the generic Least Square solver
alpha_est = utils.lsSoln(y_a,J_a,C_cal,alpha_init,epsilon,max_num_iterations,force_full_calc,plot_progress);
alpha_aoa_est = alpha_est(alpha_a_ind);
alpha_tdoa_est = alpha_est(alpha_t_ind);
alpha_fdoa_est = alpha_est(alpha_f_ind);

%% Estimate Sensor Position Error
%  Use estimated alpha, and known cal positions
%  Estimate sensor positions
beta_a_ind = 1:n_dim*n_aoa;
beta_t_ind = n_dim*n_aoa + (1:n_dim*n_tdoa);
beta_fx_ind = n_dim*(n_aoa+n_tdoa) + (1:n_dim*n_fdoa);
beta_fv_ind = n_dim*(n_aoa+n_tdoa+n_fdoa) + (1:n_dim*n_fdoa);
y_b = @(beta) reshape(z_cal - hybrid.measurement(reshape(beta(beta_a_ind),n_dim,n_aoa), ...   % x_aoa
                                                 reshape(beta(beta_t_ind),n_dim,n_tdoa), ...  % x_tdoa
                                                 reshape(beta(beta_fx_ind),n_dim,n_fdoa), ... % x_fdoa
                                                 reshape(beta(beta_fv_ind),n_dim,n_fdoa), ... % v_fdoa
                                                 x_cal, ...            % x
                                                 tdoa_ref_idx, fdoa_ref_idx, do2DAoA, ...
                                                 alpha_aoa_est, alpha_tdoa_est, alpha_fdoa_est), ...
                                                 (m_aoa+m_tdoa+m_fdoa)*n_cal,1);

J_b = @(beta) reshape(hybrid.grad_b(reshape(beta(beta_a_ind),n_dim,n_aoa), ...   % x_aoa
                                    reshape(beta(beta_t_ind),n_dim,n_tdoa), ...  % x_tdoa
                                    reshape(beta(beta_fx_ind),n_dim,n_fdoa), ... % x_fdoa
                                    reshape(beta(beta_fv_ind),n_dim,n_fdoa), ... % v_fdoa
                                    x_cal, ...            % x
                                    tdoa_ref_idx, fdoa_ref_idx, do2DAoA, ...
                                    alpha_aoa_est, alpha_tdoa_est, alpha_fdoa_est), ...
                                    n_dim*(n_aoa+n_tdoa+2*n_fdoa),(m_aoa+m_tdoa+m_fdoa)*n_cal);

% Build the initial vector
beta_init = [x_aoa(:); x_tdoa(:); x_fdoa(:); v_fdoa(:)];

beta_est = utils.lsSoln(y_b,J_b,C_cal,beta_init,epsilon,max_num_iterations,force_full_calc,plot_progress);
x_aoa_est = reshape(beta_est(beta_a_ind),n_dim,n_aoa);
x_tdoa_est = reshape(beta_est(beta_t_ind),n_dim,n_tdoa);
x_fdoa_est = reshape(beta_est(beta_fx_ind),n_dim,n_fdoa);
v_fdoa_est = reshape(beta_est(beta_fv_ind),n_dim,n_fdoa);

%% Estimate Target Position
%  Use alpha_est and beta_est with tgt measurements z to estimate

y = @(x) z - hybrid.measurement(x_aoa_est, x_tdoa_est, x_fdoa_est, v_fdoa_est,...
                                x,tdoa_ref_idx, fdoa_ref_idx, do2DAoA, ...
                                alpha_aoa_est, alpha_tdoa_est, alpha_fdoa_est);
J = @(x) hybrid.grad_x(x_aoa_est, x_tdoa_est, x_fdoa_est, v_fdoa_est,...
                       x, tdoa_ref_idx, fdoa_ref_idx, do2DAoA, ...
                       alpha_aoa_est, alpha_tdoa_est, alpha_fdoa_est);

[x,x_full] = utils.lsSoln(y,J,C_tilde,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress);