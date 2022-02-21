function [x,x_full,alpha_est,beta_est] = lsSolnUnc(x_fdoa, v_fdoa, z,C,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress,fdoa_ref_idx)
% [x,x_full,alpha_est,beta_est] = lsSolnUnc(x_fdoa, v_fdoa, z,C,x_init,epsilon,max_num_iterations,force_full_calc,
%                        plot_progress,fdoa_ref_idx)
%
% Computes the least square solution for hybrid FDOA processing.
%
% Inputs:   
%   x_fdoa              FDOA sensor positions
%   v_fdoa              FDOA sensor velocities
%   z                   Measurement vector
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
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA
%
% Outputs:
%   x               Estimated source position
%   x_full          Iteration-by-iteration estimated source positions
%   alpha_est       Array with bias estimates
%   beta_est        Array with estimated sensor positions
%
% Nicholas O'Donoughue
% 1 Nov 2021

% Parse inputs
if nargin < 10 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

% Parse inputs sizes
n_dim = size(x_fdoa,1);
n_fdoa = size(x_fdoa,2);

[test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(fdoa_ref_idx,n_fdoa);

% Initialize measurement error and Jacobian function handles
% theta vector contains x, alpha, and beta.  Let's define the
% indices
x_ind = 1:n_dim;
alpha_ind = x_ind(end) + (1:n_fdoa);
beta_fx_ind = alpha_ind(end) + (1:n_dim*n_fdoa);
beta_fv_ind = beta_fx_ind(end)+ (1:n_dim*n_fdoa);
beta_ind = [beta_fx_ind,beta_fv_ind];

y = @(theta) z - fdoa.measurement(reshape(theta(beta_fx_ind),n_dim,n_fdoa), ... % x_fdoa
                                  reshape(theta(beta_fv_ind),n_dim,n_fdoa), ... % v_fdoa
                                  reshape(theta(x_ind),n_dim,1), ...            % x
                                  fdoa_ref_idx, ...
                                  reshape(theta(alpha_ind),n_fdoa,1));   % alpha_fdoa
J = @(theta) fdoa.jacobianUnc(reshape(theta(beta_fx_ind),n_dim,n_fdoa), ... % x_fdoa
                              reshape(theta(beta_fv_ind),n_dim,n_fdoa), ... % v_fdoa
                              reshape(theta(x_ind),n_dim,1), ...            % x
                              fdoa_ref_idx, ...
                              reshape(theta(alpha_ind),n_fdoa,1));   % alpha_fdoa

% For now, we assume the AOA is independent of TDOA/FDOA
C_tilde = utils.resampleCovMtx(C, test_idx_vec, ref_idx_vec);

% Build the initial theta vector
th_init = [x_init; zeros(n_fdoa,1); x_fdoa(:); v_fdoa(:)];

% Call the generic Least Square solver
[th,th_full] = utils.lsSoln(y,J,C_tilde,th_init,epsilon,max_num_iterations,force_full_calc,plot_progress);

% Grab the x coordinates
x = th(x_ind);
x_full = th_full(x_ind,:);
alpha_est = th(alpha_ind);
beta_est = th(beta_ind);
