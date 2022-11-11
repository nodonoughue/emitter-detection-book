function [x_est,alpha_est, beta_est, A,x_grid] = mlSolnCal(x_aoa,zeta,x_cal,zeta_cal,C,C_beta,th_ctr,search_size,epsilon)
% [x_est, alpha_est, beta_est, A,x_grid] = mlSolnUnc(x_aoa,zeta,C,...
%                   x_ctr,search_size,epsilon)
%
% Construct the ML Estimate under measurement bias and sensor position
% uncertainty, by using calibration emitters and measurements to first
% estimate sensor measurement biases and positions.
%
% INPUTS:
%   x_aoa           AOA sensor positions [m]
%   zeta            Combined measurement vector
%   x_cal           Calibration sensor positions (n_dim x n_cal)
%   zeta_cal        Calibration measurements (n_aoa x n_cal)
%   C               Combined measurement error covariance matrix
%   x_ctr           Center of search grid [m]
%   search_size     Vector of search grid sizes [m]
%   epsilon         Desired resolution of search grid [m]
%
% OUTPUTS:
%   x_est           Estimated source position [m]
%   alpha_est       Array with bias estimates
%   beta_est        Array with estimated sensor positions
%   A               Likelihood computed across the entire set of candidate
%                   source positions
%   x_grid          Candidate source positions
%
% Nicholas O'Donoughue
% 16 Feb 2022


% Parse inputs sizes
n_dim = size(x_aoa,1);
n_aoa = size(x_aoa,2);

assert(size(C,1) == n_aoa || size(C,1) == 2*n_aoa,'Unable to determine if AOA measurements are 1D or 2D');
do2DAoA = size(C,1) == 2*n_aoa;
if do2DAoA
    m_aoa = 2*n_aoa;
else
    m_aoa = n_aoa;
end

%% Parse Search Space
n_th = n_dim + m_aoa + n_dim * n_aoa;
x_ind = 1:n_dim;
a_ind = n_dim + 1:m_aoa;
b_ind = a_ind(end) + 1: n_dim * n_aoa;

if isempty(th_ctr) || isscalar(th_ctr)
    % x_ctr is empty, use the origin
    x_ctr = zeros(n_dim,1);
    alpha_ctr = zeros(m_aoa,1);
    beta_ctr = x_aoa(:);

elseif numel(th_ctr) == n_dim
    % x_ctr is just the center of the target position
    
    alpha_ctr = zeros(m_aoa,1); % alpha
    beta_ctr = x_aoa(:); % beta
elseif numel(th_ctr) == n_th
    % the user already specified x_ctr to include alpha and beta
    x_ctr = th_ctr(x_ind);
    alpha_ctr = th_ctr(a_ind);
    beta_ctr = th_ctr(b_ind);
else
    error('Unable to parse center of search grid, unexpected number of dimensions received.');
end

if isempty(search_size) || isscalar(search_size)
    % search_size is empty, use default of
    % 10 km for target dimensions
    % 5 deg for angle bias
    % 25 m for sensor positions

    if isscalar(search_size)
        x_srch = search_size*ones(n_dim,1);
    else
        x_srch = 10e3*ones(n_dim,1);
    end
    alpha_srch = 5*ones(m_aoa,1);
    beta_srch = 25*ones(numel(x_aoa),1);

elseif numel(th_ctr) == n_dim
    % x_ctr is just the center of the target position
    
    x_srch = search_size(x_ind);
    alpha_srch = 5*ones(m_aoa,1);
    beta_srch = 25*ones(numel(x_aoa),1);
elseif numel(th_ctr) == n_th
    % the user already specified x_ctr to include alpha and beta
    x_srch = search_size(x_ind);
    alpha_srch = search_size(a_ind);
    beta_srch = search_size(b_ind);
else
    error('Unable to parse center of search grid, unexpected number of dimensions received.');
end

if isempty(epsilon) || isscalar(epsilon)
    % epsilon is empty, compute it based on a default desire of 101 sample
    % points per dimension for target position, and 51 for bias and sensor
    % position/velocity.
    if isscalar(epsilon)
        eps_x = epsilon * ones(n_dim,1);
    else
        eps_x = 2 * x_srch / 100;
    end
    eps_a = 2 * alpha_srch / 50;
    eps_b = 2 * beta_srch / 50;

elseif numel(epsilon) == n_dim
    % epsilon is manually defined for target position

    % computer number of elements for each dimension.
    % default is 51 per dimension for bias and sensor position/velocity.
    eps_x = epsilon;
    eps_a = 2 * alpha_srch / 50;
    eps_b = 2 * alpha_srch / 50;
elseif numel(epsilon) == n_th
    % the user already specified epsilon to include alpha and beta,
    % do nothing.
    eps_x = epsilon(x_ind);
    eps_a = epsilon(a_ind);
    eps_b = epsilon(b_ind);
else
    error('Unable to parse center of search grid, unexpected number of dimensions received.');
end

%% Estimate Measurement Bias
%  Assume target position is known (xt=x_cal(:))
%  Assume sensor position is true (beta=x_aoa)
%  Estimate unknown sensor biases (alpha_aoa)
theta_a = @(alpha) cat(1,x_cal(:),alpha,x_aoa(:));
ell_a = @(alpha) triang.loglikelihoodUnc(x_aoa, zeta_cal, C, C_beta, theta_a(alpha));

alpha_est = utils.mlSoln(ell_a,alpha_ctr,alpha_srch,eps_a);

%% Estimate Source Positions
%  Assume target position is known (xt=x_cal(:))
%  Assume sensor bias is known (previously estimated alpha_est)
%  Estimate unknown sensor positions beta
theta_b = @(beta) cat(1,x_cal(:),alpha_est,beta);
ell_b = @(beta) triang.loglikelihoodUnc(x_aoa, zeta_cal, C, C_beta, theta_b(beta));

beta_est = utils.mlSoln(ell_b,beta_ctr,beta_srch,eps_b);

%% Estimate Target Positions
%  This time, we assume that we know both alpha and beta, so we only need
%  to find tgt position x
theta_x = @(x) cat(1,x,alpha_est,beta_est);
ell = @(x) triang.loglikelihoodUnc(x_aoa, zeta, C, C_beta, theta_x(x));

[x_est,A,x_grid] = utils.mlSoln(ell,x_ctr,x_srch,eps_x);
