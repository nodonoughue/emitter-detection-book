function [x,x_full,alpha_est,beta_est] = lsSolnUnc(x_tdoa, z,C,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress,tdoa_ref_idx)
% [x,x_full,alpha_est,beta_est] = lsSolnUnc(x_tdoa, z,C,x_init,epsilon,max_num_iterations,
%                   force_full_calc, plot_progress,tdoa_ref_idx)
%
% Computes the gradient descent solution for TDOA processing.
%
% Inputs:   
%   x_tdoa              TDOA sensor positions
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
%   tdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for TDOA
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
if nargin < 11 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

% Parse inputs sizes
n_dim = size(x_tdoa,1);
n_tdoa = size(x_tdoa,2);

[test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(tdoa_ref_idx,n_tdoa);

% Initialize measurement error and Jacobian function handles
% theta vector contains x, alpha, and beta.  Let's define the
% indices
x_ind = 1:n_dim;
alpha_ind = x_ind(end) + (1:n_tdoa);
beta_ind = alpha_ind(end) + (1:n_dim*n_tdoa);

y = @(theta) z - tdoa.measurement(reshape(theta(beta_ind),n_dim,n_tdoa), ...  % x_tdoa
                                  reshape(theta(x_ind),n_dim,1), ...            % x
                                  tdoa_ref_idx, ...
                                  reshape(theta(alpha_ind),n_tdoa,1));   % alpha_tdoa
J = @(theta) tdoa.jacobianUnc(reshape(theta(beta_ind),n_dim,n_tdoa), ...  % x_tdoa
                              reshape(theta(x_ind),n_dim,1), ...            % x
                              tdoa_ref_idx, ...
                              reshape(theta(alpha_ind),n_tdoa,1)); % alpha_tdoa

% Parse the TDOA and FDOA reference indices together
C_tilde = utils.resampleCovMtx(C, test_idx_vec, ref_idx_vec);

% Build the initial theta vector
th_init = [x_init; zeros(n_tdoa,1); x_tdoa(:)];

% Call the generic Least Square solver
[th,th_full] = utils.lsSoln(y,J,C_tilde,th_init,epsilon,max_num_iterations,force_full_calc,plot_progress);

% Grab the x coordinates
x = th(x_ind);
x_full = th_full(x_ind,:);
alpha_est = th(alpha_ind);
beta_est = th(beta_ind);

