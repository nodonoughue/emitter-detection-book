function ell = loglikelihood(x_fdoa,v_fdoa,rho_dot,C,x_source,ref_idx)
% ell = loglikelihood(x_fdoa,v_fdoa,rho_dot,C,x_source,ref_idx)
%
% Computes the Log Likelihood for FDOA sensor measurement, given the 
% received measurement vector rho_dot, covariance matrix C, 
% and set of candidate source positions x_source.
%
% INPUTS:
%   x_fdoa      Sensor positions [m]
%   v_fdoa      Sensor velocities [m/s]
%   rho_dot     FDOA measurement vector
%   C           FDOA measurement error covariance matrix
%   x_source    Candidate source positions
%   ref_idx     Scalar index of reference sensor, or nDim x nPair
%               matrix of sensor pairings
%
% OUTPUTS:
%   ell         Log-likelihood evaluated at each position x_source.
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 6 || ~exist('ref_idx','var')
    ref_idx = [];
end

n_source_pos = size(x_source,2);
ell = zeros(n_source_pos,1);

% Preprocess covariance matrix
C_d = decomposition(C);

for idx_source = 1:n_source_pos
    % Generate the ideal measurement matrix for this position
    r_dot = fdoa.measurement(x_fdoa,v_fdoa, x_source(:,idx_source), ref_idx);
    
    % Compute the log-likelihood
    err = (rho_dot - r_dot);
    
    ell(idx_source) = -err'/C_d*err;
end