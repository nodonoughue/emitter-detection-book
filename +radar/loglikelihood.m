function ell = loglikelihood(x_rdr,v_rdr,zeta,C,x_tgt,v_tgt,angle_dims)
% ell = loglikelihood(x_rdr,v_rdr,zeta,C,x_tgt,v_tgt,angle_dims)
%
% Computes the Log Likelihood for monostatic range measurements, given
% the received measurement vector rho, covariance matrix C, 
% and set of candidate target positions x_tgt and velocities v_tgt.
%
% The target positions and velocities must be broadcastable.  To have
% paired pos/vel entries, supply both as nDim x nTgts, but to test all
% possible pairs, supply them as nDim x nPos x 1 (position) and
% nDim x 1 x nVel (velocities).
%
% Inputs:
%   
%   x_rdr               Radar positions [m]
%   v_rdr               Radar velocities [m/s]
%   zeta                Measurements [m, m/s, rad]
%   C                   Measurement Error Covariance Matrix [m^2, m^2/s^2,
%                       rad^2]
%   x_tgt               Candidate target position [m]
%   v_tgt               Candidate target velocity [m/s]
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
% 15 September 2022

%% Parse inputs
if nargin < 7 || ~exist('angle_dims','var') || isempty(angle_dims)
    angle_dims = 0;
end
assert(angle_dims >=0 && angle_dims <= 2,'Error parsing angle dimensions command, must be between 0 and 2.');

do_doppler = ~isempty(v_rdr) || ~isempty(v_rx) || ~isempty(v_tgt);

if do_doppler && isempty(v_rdr)
    v_rdr = zeros(size(x_rdr));
end

if do_doppler && isempty(v_tgt)
    v_tgt = zeros(size(x_tgt));
end

% Check source pos/vel inputs
if do_doppler
    size_pos = size(x_tgt);
    size_vel = size(v_tgt);
    assert( ~any((size_pos ~= size_vel) & (size_pos > 1) & (size_vel > 1)), 'Non-singleton dimensions of source position and velocity must match.');

    in_size = max(size_pos,size_vel);
    out_size = in_size(2:end); % drop first dimension
    nDim = in_size(1);
    n_source_pos = prod(out_size);

    % Broadcast x_source and v_source to common size
    x_tgt = reshape(x_tgt.*ones(size(in_size)), nDim, n_source_pos);
    v_tgt = reshape(v_tgt.*ones(size(in_size)), nDim, n_source_pos);
else
    in_size = size(x_tgt);
    out_size = in_size(2:end); % drop first dimension
    nDim = in_size(1);
    n_source_pos = prod(out_size);

    x_tgt = reshape(x_tgt,nDim,n_source_pos);
    v_tgt = zeros(nDim,n_source_pos);
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
    x_i = x_tgt(:,idx_source);
    v_i = v_tgt(:,idx_source);
    
    % Generate the ideal measurement matrix for this position
    if do_doppler
        r = radar.measurement(x_rdr, x_i, v_rdr, v_i, angle_dims);
    else
        r = radar.measurement(x_rdr, x_i, [], [], angle_dims);
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