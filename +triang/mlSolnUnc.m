function [x_est,alpha_est, beta_est, A,x_grid] = mlSolnUnc(x_aoa,zeta,C,x_ctr,search_size,epsilon)
% [x_est, alpha_est, beta_est, A,x_grid] = mlSolnUnc(x_aoa,zeta,C,...
%                   x_ctr,search_size,epsilon)
%
% Construct the ML Estimate by systematically evaluating the log
% likelihood function at a series of coordinates, and returning the index
% of the maximum.  Optionally returns the full set of evaluated
% coordinates, as well.
%
% INPUTS:
%   x_aoa           AOA sensor positions [m]
%   zeta            Combined measurement vector
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

% Initialize measurement error and Jacobian function handles
% theta vector contains x, alpha, and beta.  Let's define the
% indices
x_ind = 1:n_dim;
a_ind = n_dim + (1:m_aoa);
b_ind = a_ind(end) + 1:(n_dim*n_aoa);

% Set up function handle
% We must take care to ensure that it can handle an n_th x N matrix of
% inputs; for compatibility with how utils.mlSoln will call it.
ell = @(theta) triang.loglikelihoodUnc(x_aoa, zeta, C, theta);
        
%% Define center of search space
n_th = n_dim + m_aoa + n_dim * n_aoa;
if isempty(x_ctr) || isscalar(x_ctr)
    % x_ctr is empty, use the origin

    th_ctr = [zeros(n_dim,1);
              zeros(m_aoa,1); % alpha
              x_aoa(:)]; % beta
elseif numel(x_ctr) == n_dim
    % x_ctr is just the center of the target position
    
    % build th_ctr manually
    th_ctr = [x_ctr; 
              zeros(m_aoa,1); % alpha
              x_aoa(:)]; % beta
elseif numel(x_ctr) == n_th
    % the user already specified x_ctr to include alpha and beta
    th_ctr = x_ctr;
else
    error('Unable to parse center of search grid, unexpected number of dimensions received.');
end

%% Define search size
if isempty(search_size) || isscalar(search_size)
    % search_size is empty, use default of
    % 10 km for target dimensions
    % 1 deg for angle bias
    % 25 m for sensor positions

    search_size = [10e3 * ones(1, n_dim), ...      % x [m]
                   pi/180 * ones(1, m_aoa), ...    % alpha_aoa [rad]
                   25 * ones(1, n_aoa*n_dim)];     % x_aoa [m]
elseif numel(search_size) == n_dim
    % search_size is just defined for the target position
    
    % build th_ctr manually
    search_size = [search_size(:)', ...            % x [m]
                   pi/180 * ones(1, m_aoa), ...    % alpha_aoa [rad]
                   25 * ones(1, n_aoa*n_dim)];     % x_aoa [m]
elseif numel(search_size) == n_th
    % the user already specified search_size to include alpha and beta,
    % do nothing.
else
    error('Unable to parse center of search grid, unexpected number of dimensions received.');
end

%% Define grid resolution
if isempty(epsilon) || isscalar(epsilon)
    % epsilon is empty, compute it based on a default desire of 101 sample
    % points per dimension for target position, and 11 for bias and sensor
    % position/velocity.
    n_elements = [101 * ones(1, n_dim), ...
                  11 * ones(1, m_aoa), ...
                  11 * ones(1, n_dim * n_aoa )];

    epsilon = 2*search_size ./ (n_elements-1);

elseif numel(epsilon) == n_dim
    % epsilon is manually defined for target position

    % computer number of elements for each dimension.
    % default is 11 per dimension for bias and sensor position/velocity.
    n_elements = [1 + 2*search_size(x_ind)./epsilon, ...
                  11 * ones(1, m_aoa), ...
                  11 * ones(1, n_dim * n_aoa)];
   
    % re-compute grid spacing
    epsilon = 2*search_size ./ (n_elements-1);

elseif numel(epsilon) == n_th
    % the user already specified epsilon to include alpha and beta,
    % do nothing.
else
    error('Unable to parse center of search grid, unexpected number of dimensions received.');
end

%% Call the util function
[th_est,A,x_grid] = utils.mlSoln(ell,th_ctr,search_size,epsilon);

x_est = th_est(x_ind);
alpha_est = th_est(a_ind);
beta_est = th_est(b_ind);