function [x_est,A,x_grid] = bfSolnUnc(x_aoa,x_tdoa,x_fdoa,v_fdoa,zeta,C,x_ctr,search_size,epsilon,tdoa_ref_idx,fdoa_ref_idx,pdftype)
% [x_est,A,x_grid] = bfSolnUnc(x_aoa,x_tdoa,x_fdoa,v_fdoa,zeta,C,x_ctr,
%                              search_size,epsilon,tdoa_ref_idx,
%                              fdoa_ref_idx,pdftype)
%
% Construct the BestFix estimate by systematically evaluating the PDF at 
% a series of coordinates, and returning the index of the maximum.  
% Optionally returns the full set of evaluated coordinates, as well.
%
% Assumes a multi-variate Gaussian distribution with covariance matrix C,
% and unbiased estimates at each sensor.  Note that the BestFix algorithm
% implicitly assumes each measurement is independent, so any cross-terms in
% the covariance matrix C are ignored.
%
% INPUTS:
%   x_aoa           AOA sensor positions [m]
%   x_tdoa          TDOA sensor positions [m]
%   x_fdoa          FDOA sensor positions [m]
%   v_fdoa          FDOA sensor velocities [m/s]
%   zeta            Combined measurement vector
%   C               Combined measurement error covariance matrix
%   x_ctr           Center of search grid [m]
%   search_size     2-D vector of search grid sizes [m]
%   epsilon         Desired resolution of search grid [m]
%   tdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for TDOA
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA
%   pdftype         String indicating the type of distribution to use.
%                   See +utils/makePDFs.m for options.
%
% OUTPUTS:
%   x_est           Estimated source position [m]
%   A               Likelihood computed across the entire set of candidate
%                   source positions
%   x_grid          Candidate source positions
%
% Ref:
%  Eric Hodson, "Method and arrangement for probabilistic determination of 
%  a target location," U.S. Patent US5045860A, 1990, 
%  https://patents.google.com/patent/US5045860A
%
%
% Nicholas A. O'Donoughue
% 1 Nov 2021

if nargin < 10 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

if nargin < 11 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

if nargin < 12 || isempty(pdftype)
    pdftype = [];
end

% Parse inputs sizes
n_dim = size(x_aoa,1);
n_aoa = size(x_aoa,2);
n_tdoa = size(x_tdoa,2);
n_fdoa = size(x_fdoa,2);


assert(size(C,1) == n_aoa + n_tdoa + n_fdoa || size(C,1) == 2*n_aoa + n_tdoa + n_fdoa,'Unable to determine if AOA measurements are 1D or 2D');
do2DAoA = size(C,1) == 2*n_aoa + n_tdoa + n_fdoa;
if do2DAoA
    m_aoa = 2*n_aoa;
else
    m_aoa = n_aoa;
end

[tdoa_test_idx_vec, tdoa_ref_idx_vec] = utils.parseReferenceSensor(tdoa_ref_idx,n_tdoa);
[fdoa_test_idx_vec, fdoa_ref_idx_vec] = utils.parseReferenceSensor(fdoa_ref_idx,n_fdoa);

m_tdoa = numel(tdoa_test_idx_vec);
m_fdoa = numel(fdoa_test_idx_vec);

% Initialize measurement error and Jacobian function handles
% theta vector contains x, alpha, and beta.  Let's define the
% indices
x_ind = 1:n_dim;
alpha_a_ind = x_ind(end) + (1:m_aoa);
alpha_t_ind = alpha_a_ind(end) + (1:m_tdoa);
alpha_f_ind = alpha_t_ind(end) + (1:m_fdoa);
beta_a_ind = alpha_f_ind(end) + (1:n_dim*n_aoa);
beta_t_ind = beta_a_ind(end) + (1:n_dim*n_tdoa);
beta_fx_ind = beta_t_ind(end) + (1:n_dim*n_fdoa);
beta_fv_ind = beta_fx_ind(end)+ (1:n_dim*n_fdoa);

% Parse the TDOA and FDOA reference indices together
test_idx_vec = cat(2,tdoa_test_idx_vec, n_tdoa + fdoa_test_idx_vec);
ref_idx_vec = cat(2,tdoa_ref_idx_vec, n_tdoa + fdoa_ref_idx_vec);

% For now, we assume the AOA is independent of TDOA/FDOA
C_aoa = C(1:m_aoa, 1:m_aoa);
C_tfdoa = C(m_aoa+1:end, m_aoa+1:end);
C_tilde = blkdiag(C_aoa, utils.resampleCovMtx(C_tfdoa, test_idx_vec, ref_idx_vec));

% Generate the PDF
z = @(theta) hybrid.measurement(reshape(theta(beta_a_ind),n_dim,n_aoa), ...   % x_aoa
                                reshape(theta(beta_t_ind),n_dim,n_tdoa), ...  % x_tdoa
                                reshape(theta(beta_fx_ind),n_dim,n_fdoa), ... % x_fdoa
                                reshape(theta(beta_fv_ind),n_dim,n_fdoa), ... % v_fdoa
                                reshape(theta(x_ind),n_dim,1), ...            % x
                                tdoa_ref_idx, fdoa_ref_idx, ...
                                reshape(theta(alpha_a_ind),n_aoa,1), ...      % alpha_aoa
                                reshape(theta(alpha_t_ind),m_tdoa,1), ...% alpha_tdoa
                                reshape(theta(alpha_f_ind),m_fdoa,1));   % alpha_fdoa

pdfs = utils.makePDFs(z, zeta, pdftype, C_tilde);


%% Define center of search space
n_th = n_dim + m_aoa + m_tdoa + m_fdoa + n_dim * (n_aoa + n_tdoa + 2*n_fdoa);
if isempty(x_ctr) || isscalar(x_ctr)
    % x_ctr is empty, use the origin

    th_ctr = [zeros(n_dim,1);
              zeros(m_aoa + m_tdoa + m_fdoa,1); % alpha
              x_aoa(:); x_tdoa(:); x_fdoa(:); v_fdoa(:)]; % beta
elseif numel(x_ctr) == n_dim
    % x_ctr is just the center of the target position
    
    % build th_ctr manually
    th_ctr = [x_ctr; 
              zeros(m_aoa + m_tdoa + m_fdoa,1); % alpha
              x_aoa(:); x_tdoa(:); x_fdoa(:); v_fdoa(:)]; % beta
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
    % 10 m for range bias
    % 10 m/s for range-rate bias
    % 25 m for sensor positions
    % 10 m/s for sensor velocities

    search_size = [10e3 * ones(1, n_dim), ...      % x [m]
                   pi/180 * ones(1,m_aoa), ...   % alpha_aoa [rad]
                   10 * ones(1,m_tdoa), ...        % alpha_tdoa [m]
                   10 * ones(1,m_fdoa), ...        % alpha_fdoa [m/s]
                   25 * ones(1, n_aoa*n_dim), ...  % x_aoa [m]
                   25 * ones(1, n_tdoa*n_dim), ... % x_tdoa [m]
                   25 * ones(1, n_fdoa*n_dim), ... % x_fdoa [m]
                   10 * ones(1, n_fdoa*n_dim)];    % v_fdoa [m/s]           
elseif numel(search_size) == n_dim
    % search_size is just defined for the target position
    
    % build th_ctr manually
    search_size = [search_size(:)', ...            % x [m]
                   pi/180 * ones(1,m_aoa), ...   % alpha_aoa [rad]
                   10 * ones(1,m_tdoa), ...        % alpha_tdoa [m]
                   10 * ones(1,m_fdoa), ...        % alpha_fdoa [m/s]
                   25 * ones(1, n_aoa*n_dim), ...  % x_aoa [m]
                   25 * ones(1, n_tdoa*n_dim), ... % x_tdoa [m]
                   25 * ones(1, n_fdoa*n_dim), ... % x_fdoa [m]
                   10 * ones(1, n_fdoa*n_dim)];    % v_fdoa [m/s]      
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
                  11 * ones(1, m_aoa + m_tdoa + m_fdoa), ...
                  11 * ones(1, n_dim * (n_aoa + n_tdoa + 2*n_fdoa))];

    epsilon = 2*search_size ./ (n_elements-1);

elseif numel(epsilon) == n_dim
    % epsilon is manually defined for target position

    % computer number of elements for each dimension.
    % default is 11 per dimension for bias and sensor position/velocity.
    n_elements = [1 + 2*search_size(x_ind)./epsilon, ...
                  11 * ones(1, m_aoa + m_tdoa + m_fdoa), ...
                  11 * ones(1, n_dim * (n_aoa + n_tdoa + 2*n_fdoa))];
   
    % re-compute grid spacing
    epsilon = 2*search_size ./ (n_elements-1);

elseif numel(epsilon) == n_th
    % the user already specified epsilon to include alpha and beta,
    % do nothing.
else
    error('Unable to parse center of search grid, unexpected number of dimensions received.');
end

%% Call the util function
[theta_est,A,x_grid] = utils.bestfix(pdfs,th_ctr,search_size,epsilon);

x_est = theta_est(idx_dim);