function crlb = computeCRLB(x_tdoa,xs,C,ref_idx,variance_is_toa,resample_covariance)
% crlb = computeCRLB(x_tdoa,xs,C,ref_idx,variance_is_toa,resample_covariance)
%
% Computes the CRLB on position accuracy for source at location xs and
% sensors at locations in x1 (Ndim x N).  Ctdoa is an Nx1 vector of TOA
% variances at each of the N sensors.
%
% Inputs:
%   x_tdoa      (Ndim x N) array of TDOA sensor positions
%   xs          (Ndim x M) array of source positions over which to 
%               calculate CRLB
%   C           TOA covariance matrix [s^2]
%   ref_idx     Scalar index of reference sensor, or nDim x nPair
%               matrix of sensor pairings
%   variance_is_toa (Optional) flag indicating whether supplied variance is
%               in units of time (TRUE) or distance (FALSE). Default=TRUE.
%   resample_covariance (Optional) flag indicating whether the covariance 
%               matrix should be resamples (TRUE) to convert from sensor
%               errors to measurement errors, or whether it was directly
%               supplied as measurement errors (FALSE). Default=TRUE.
%
% Outputs:
%   crlb    Lower bound on the error covariance matrix for an unbiased
%           TDOA estimator (Ndim x Ndim)
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 6 || ~exist('resample_covariance','var')
    resample_covariance = true;
end

if nargin < 5 || ~exist('variance_is_toa','var')
    variance_is_toa = true;
end

if nargin < 4 || ~exist('ref_idx','var')
    ref_idx = [];
end

[~, n_sensor] = size(x_tdoa);

% Construct Jacobian function handle
J = @(x) tdoa.jacobian(x_tdoa,x,ref_idx);

% Preprocess covariance matrix
if variance_is_toa
    C_out = C*utils.constants.c^2;
else
    C_out = C;
end

% Parse sensor pairs
if resample_covariance
    [test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(ref_idx, n_sensor);
    C_tilde = utils.resampleCovMtx(C_out, test_idx_vec, ref_idx_vec);
else
    C_tilde = C_out;
end

% Call generic solver
crlb = utils.computeCRLB(xs,C_tilde,J);