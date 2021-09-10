function ell = loglikelihood(x_sensor, psi,C,x_source,do2DAoA)
% function ell = loglikelihood(x_sensor,psi,C,x_source)
%
% Computes the Log Likelihood for AOA sensor measurement, given the 
% received measurement vector psi, covariance matrix C, 
% and set of candidate source positions x_source.
%
% If 2D AOA (az/el) is desired, then the measurement vector psi should be
% stacked (psi = [az(:); el(:)]).  In this case, the provided covariance
% matrix should have dimension 2*nSource x 2*nSource.
%
% INPUTS:
%   x_sensor    nDim x nAOA vector of AOA sensor positions
%   psi         Measurement vector
%   C           Measurement error covariance matrix
%   x_source    Candidate source positions
%   do2DAoA     Boolean flag indicating whether 2D AOA (az/el) should be
%               conducted
%
% OUTPUTS:
%   ell         Log-likelihood evaluated at each position x_source.
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
n_source_pos = size(x_source,2);
ell = zeros(n_source_pos,1);

if nargin < 5 || isempty(do2DAoA)
    do2DAoA = true;
end

if numel(psi) ~= size(C, 1) || numel(psi) ~= size(C, 2)
    error('Covariance matrix and measurement vector sizes do not match.');
end

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
    % Generate the ideal measurement matrix for this position
    p = triang.measurement(x_sensor, x_source(:,idx_source), do2DAoA);
    
    % Compute the measurement error
    err = (psi - p);

    % Evaluate the scaled log likelihood
    if do_decomp
        ell(idx_source) = -err'/C_d*err;
    else
        ell(idx_source) = -err*C_inv*err;
    end
end