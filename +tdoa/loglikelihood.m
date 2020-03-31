function ell = loglikelihood(x_sensor,rho,C,x_source,ref_idx)
% ell = loglikelihood(x_sensor,rho,C,x_source,ref_idx)
%
% Computes the Log Likelihood for TDOA sensor measurement, given the 
% received measurement vector rho, covariance matrix C, 
% and set of candidate source positions x_source.
%
% INPUTS:
%   x_sensor    TDOA sensor positions 
%   rho         Measurement vector
%   C           measurement error covariance matrix
%   x_source    Candidate source position
%   ref_idx     Scalar index of reference sensor, or nDim x nPair
%               matrix of sensor pairings
%
% OUTPUTS:
%   ell         Log-likelihood evaluated at each position x_source.
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 5 || ~exist('ref_idx','var')
    ref_idx = [];
end

n_source_pos = size(x_source,2);
ell = zeros(n_source_pos,1);

% Preprocess covariance matrix
C_d = decomposition(C);

for idx_source = 1:n_source_pos
    % Generate the ideal measurement matrix for this position
    r = tdoa.measurement(x_sensor, x_source(:,idx_source), ref_idx);
    
    % Compute the log-likelihood
    err = (rho - r);
    
    ell(idx_source) = -err'/C_d*err;
end