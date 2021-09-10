function ell = loglikelihood(x_tx,x_rx,rho,C,x_source,ref_idx)
% ell = loglikelihood(x_tx,x_rx,rho,C,x_source,ref_idx)
%
% Computes the Log Likelihood for PCL bistatic range measurements, given
% the received measurement vector rho, covariance matrix C, 
% and set of candidate source positions x_source.
%
% TODO: Incorporate range rate measurements
% TODO: Incorporate angle of arrival measurements
%
% Inputs:
%   
%   x_tx                Transmitter positions [m]
%   x_rx                Receiver positions [m]
%   rho                 Bistatic Range Measurements [m]
%   C                   Measurement Error Covariance Matrix [m^2]
%   x_source            Candidate source position
%   ref_idx             Matrix of tx/rx pairing indices (in the order
%                       they're used in C).  If ignored, then all pairwise
%                       measurements are used (nTx x nRx)
%
% OUTPUTS:
%   ell         Log-likelihood evaluated at each position x_source.
%
% Nicholas O'Donoughue
% 10 September 2021

% Parse inputs
if nargin < 6 || ~exist('ref_idx','var')
    ref_idx = [];
end

n_source_pos = size(x_source,2);
ell = zeros(n_source_pos,1);

% Pre-compute covariance matrix inverses
do_decomp = ~verLessThan('MATLAB','9.3');
if do_decomp
    % Starging in R2017b, MATLAB released the DECOMPOSITION function,
    % which can decompose matrices for faster computation of left- and
    % right-division in for loops.
    C_d = decomposition(C,'chol');
else
    % If DECOMPOSITION is unavailable, let's precompute the pseudo-inverse.
    C_inv = pinv(C);
end

for idx_source = 1:n_source_pos
    x_i = x_source(:,idx_source);
    
    % Generate the ideal measurement matrix for this position
    r = pcl.measurement(x_tx,x_rx, x_i, [],[],[],ref_idx);
    
    % Compute the measurement error
    err = (rho - r);

    % Evaluate the scaled log likelihood
    if do_decomp
        ell(idx_source) = -err'/C_d*err;
    else
        ell(idx_source) = -err*C_inv*err;
    end
end