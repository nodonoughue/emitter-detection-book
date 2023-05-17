function crlb = computeCRLBpriors(x_fdoa,v_fdoa,xs,C,fim_prior,ref_idx)
% crlb = computeCRLBpriors(x_fdoa,v_fdoa,xs,C,fim_prior,ref_idx)
%
% Computes the CRLB on position accuracy for source at location xs and
% sensors at locations in x_fdoa (Ndim x N) with velocity v_fdoa.  
% C is an Nx1 vector of FOA variances at each of the N sensors, and ref_idx
% defines the reference sensor(s) used for FDOA.
%
% This version computes the CRLB in the presence of a statistical prior,
% created using the util.makePrior function, which returns not only a
% function handle for the prior, but the FIM for the log of the prior. The
% combined CRLB is computed according to (5.38).
%
% Inputs:
%   x_fdoa      (Ndim x N) array of FDOA sensor positions
%   v_fdoa      (Ndim x N) array of FDOA sensor velocities
%   xs          (Ndim x M) array of source positions over which to 
%               calculate CRLB
%   C           Covariance matrix for range rate estimates at the N
%               FDOA sensors [(m/s)^2]
%   fim_prior   Function handle that computes the Fisher Information Matrix
%               for the statistical prior
%   ref_idx     Scalar index of reference sensor, or nDim x nPair
%               matrix of sensor pairings
%   
% Outputs:
%   crlb    Lower bound on the error covariance matrix for an unbiased
%           FDOA estimator (Ndim x Ndim)
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 5 || ~exist('ref_idx','var')
    ref_idx = [];
end
[~, n_sensor] = size(x_fdoa);

% Define the Jacobian function
J = @(x) fdoa.jacobian(x_fdoa,v_fdoa,x,ref_idx);

% Resample covariance matrix
[test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(ref_idx, n_sensor);
C_tilde = utils.resampleCovMtx(C, test_idx_vec, ref_idx_vec);

% Ensure the covariance matrix is invertible
C_tilde = utils.ensureInvertible(C_tilde);

% Call the generic solver with the prior information
crlb = utils.computeCRLB(xs, C_tilde, J, fim_prior);