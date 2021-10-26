function ell = loglikelihood(x_tx,x_rx,v_tx,v_rx,zeta,C,x_src,v_src,ref_idx,angle_dims)
% ell = loglikelihood(x_tx,x_rx,v_tx,v_rx,rho,C,x_src,v_src,ref_idx,angle_dims)
%
% Computes the Log Likelihood for PCL bistatic range measurements, given
% the received measurement vector rho, covariance matrix C, 
% and set of candidate source positions x_source and velocities v_source.
%
% The source positions and velocities must be broadcastable.  To have
% paired pos/vel entries, supply both as nDim x nSources, but to test all
% possible pairs, supply them as nDim x nPos x 1 (position) and
% nDim x 1 x nVel (velocities).
%
% Inputs:
%   
%   x_tx                Transmitter positions [m]
%   x_rx                Receiver positions [m]
%   v_tx                Transmitter velocities [m/s]
%   v_rx                Receiver velocities [m/s]
%   zeta                Bistatic Measurements [m, m/s, rad]
%   C                   Measurement Error Covariance Matrix [m^2, m^2/s^2,
%                       rad^2]
%   x_src               Candidate source position [m]
%   v_src               Candidate source velocity [m/s]
%   ref_idx             Matrix of tx/rx pairing indices (in the order
%                       they're used in C).  If ignored, then all pairwise
%                       measurements are used (nTx x nRx).
%   angle_dims          Number of receiver angle of arrival dimensions
%                       reported (0 = no angles, 1 = azimuth, 2 = azimuth 
%                       & elevation)
%
% OUTPUTS:
%   ell         Log-likelihood evaluated at each position x_source and
%               velocity v_source.  Shape will be the common shape of
%               x_source and v_source (ignoring the first dimension).
%
% Nicholas O'Donoughue
% 26 October 2021

%% Parse inputs
if nargin < 10 || ~exist('angle_dims','var') || isempty(angle_dims)
    angle_dims = 0;
end
assert(angle_dims >=0 && angle_dims <= 2,'Error parsing angle dimensions command, must be between 0 and 2.');

do_doppler = ~isempty(v_tx) || ~isempty(v_rx) || ~isempty(v_src);

if nargin < 9 || ~exist('ref_idx','var')
    ref_idx = [];
end

if do_doppler && isempty(v_tx)
    v_tx = zeros(size(x_tx));
end

if do_doppler && isempty(v_rx)
    v_rx = zeros(size(x_rx));
end

if do_doppler && isempty(v_src)
    v_src = zeros(size(x_tgt));
end

% Check source pos/vel inputs
if do_doppler
    size_pos = size(x_src);
    size_vel = size(v_src);
    assert( ~any((size_pos ~= size_vel) & (size_pos > 1) & (size_vel > 1)), 'Non-singleton dimensions of source position and velocity must match.');

    in_size = max(size_pos,size_vel);
    out_size = in_size(2:end); % drop first dimension
    nDim = in_size(1);
    n_source_pos = prod(out_size);

    % Broadcast x_source and v_source to common size
    x_src = reshape(x_src.*ones(size(in_size)), nDim, n_source_pos);
    v_src = reshape(v_src.*ones(size(in_size)), nDim, n_source_pos);
else
    in_size = size(x_src);
    out_size = in_size(2:end); % drop first dimension
    nDim = in_size(1);
    n_source_pos = prod(out_size);

    x_src = reshape(x_src,nDim,n_source_pos);
    v_src = zeros(nDim,n_source_pos);
end

%% Initialize Outputs
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
    x_i = x_src(:,idx_source);
    v_i = v_src(:,idx_source);
    
    % Generate the ideal measurement matrix for this position
    if do_doppler
        r = pcl.measurement(x_tx, x_rx, x_i, v_tx, v_rx, v_i, ref_idx, angle_dims);
    else
        r = pcl.measurement(x_tx, x_rx, x_i, [], [], [], ref_idx, angle_dims);
    end

    % Compute the measurement error
    err = (zeta - r);

    % Evaluate the scaled log likelihood
    if do_decomp
        ell(idx_source) = -err'/C_d*err;
    else
        ell(idx_source) = -err*C_inv*err;
    end
end