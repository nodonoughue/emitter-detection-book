function crlb = computeCRLB(x_aoa,x_tdoa,x_fdoa,v_fdoa,xs,C,tdoa_ref_idx,fdoa_ref_idx)
% crlb = computeCRLB(x_aoa,x_tdoa,x_fdoa,v_fdoa,xs,C,tdoa_ref_idx,
%                                                    fdoa_ref_idx)
%
% Computes the CRLB on position accuracy for source at location xs and
% a combined set of AOA, TDOA, and FDOA measurements.  The covariance
% matrix C dictates the combined variances across the three measurement
% types.
%
% Inputs:
%   x_aoa           (Ndim x Na) array of AOA sensor positions
%   x_tdoa          (Ndim x Nt) array of TDOA sensor positions
%   x_fdoa          (Ndim x Nf) array of FDOA sensor positions
%   v_fdoa          (Ndim x Nf) array of FDOA sensor velocities
%   xs              (Ndim x M) array of source positions over which to 
%                   calculate CRLB
%   C               Combined covariance matrix for the AOA, TDOA and FDOA
%                   estimates (Na + Nt + Nf x Na + Nt + Nf)
%   tdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for TDOA
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA
%
% Outputs:
%   crlb    Lower bound on the error covariance matrix for an unbiased
%           TDOA estimator (Ndim x Ndim)
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 6 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end
if nargin < 7 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end
[n_dim,n_source] = size(xs);

% Initialize Jacobian
J = @(x) hybrid.jacobian(x_aoa,x_tdoa,x_fdoa,v_fdoa,x,tdoa_ref_idx,fdoa_ref_idx);

% Resample the covariance matrix
n_aoa = size(x_aoa,2);
n_tdoa = size(x_tdoa,2);
n_fdoa = size(x_fdoa,2);

% Parse the TDOA and FDOA reference indices together
[tdoa_test_idx_vec, tdoa_ref_idx_vec] = utils.parseReferenceSensor(tdoa_ref_idx,n_tdoa);
[fdoa_test_idx_vec, fdoa_ref_idx_vec] = utils.parseReferenceSensor(fdoa_ref_idx,n_fdoa);
test_idx_vec = [1:n_aoa, n_aoa + tdoa_test_idx_vec, n_aoa + n_tdoa + fdoa_test_idx_vec];
ref_idx_vec = [nan*ones(1,n_aoa), n_aoa + tdoa_ref_idx_vec, n_aoa + n_tdoa + fdoa_ref_idx_vec];

% Resample the covariance matrix
C_tilde = utils.resampleCovMtx(C, test_idx_vec, ref_idx_vec);

% Make sure it's invertible
C_tilde = utils.ensureInvertible(C_tilde);

% Pre-compute covariance matrix inverses
do_decomp = ~verLessThan('MATLAB','9.3');
if do_decomp
    % Starging in R2017b, MATLAB released the DECOMPOSITION function,
    % which can decompose matrices for faster computation of left- and
    % right-division in for loops.
    C_d = decomposition(C_tilde,'chol');
else
    % If DECOMPOSITION is unavailable, let's precompute the pseudo-inverse.
    C_inv = pinv(C_tilde);
end

% Initialize output variable
crlb = zeros([n_dim,n_dim,n_source]);

% Repeat CRLB for each of the n_source test positions
for idx =1:n_source
    this_x = xs(:,idx);
    
    % Evaluate the Jacobian
    J_i = J(this_x);
    
    % Compute the Fisher Information Matrix
    if do_decomp
        F = J_i/C_d*J_i'; % Ndim x Ndim
    else
        F = J_i*C_inv*J_i';
    end
    
    if any(isnan(F(:))) || any(isinf(F(:)))
        % Configuration is ill-posed; likely because source position
        % overlaps with the reference sensor
        crlb(:,:,idx) = NaN;
    else
        % Invert the Fisher Information Matrix to compute the CRLB
        crlb(:,:,idx) = pinv(F);
    end
end

