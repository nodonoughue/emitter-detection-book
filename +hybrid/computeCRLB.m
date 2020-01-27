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

% Preprocess covariance matrix
C_d = decomposition(C);

% Initialize output variable
crlb = zeros([n_dim,n_dim,n_source]);

% Repeat CRLB for each of the n_source test positions
for idx =1:n_source
    this_x = xs(:,idx);
    J_i = J(this_x);
    F = J_i/C_d*J_i'; % Ndim x Ndim
    if any(isnan(F(:))) || any(isinf(F(:)))
        % Configuration is ill-posed; likely because source position
        % overlaps with the reference sensor
        crlb(:,:,idx) = NaN;
    else
        % Invert the Fisher Information Matrix to compute the CRLB
        crlb(:,:,idx) = pinv(F);
    end
end

