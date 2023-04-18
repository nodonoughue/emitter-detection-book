function crlb = computeCRLB(x_tx,x_rx,x_src,v_tx,v_rx,v_src,C,ref_idx,angle_dims)
% crlb = computeCRLB(x_tx,x_rx,x_src,v_tx,v_rx,v_src,C,ref_idx,angle_dims)
%
% Computes the CRLB on position accuracy for source at location x_src and
% velocity v_src, given sensors at locations x_rx with velocities v_rx,
% and transmitters at x_tx with velocities v_tx. Assumes perfect knowledge
% of sensor and transmitter parameters.
%
% Inputs:
%   
%   x_tx                Transmitter positions [m]
%   x_rx                Receiver positions [m]
%   x_src               Candidate source position [m]
%   v_tx                Transmitter velocities [m/s]
%   v_rx                Receiver velocities [m/s]
%   v_src               Candidate source velocity [m/s]
%   C                   Measurement Error Covariance Matrix [m^2, m^2/s^2,
%                       rad^2]
%   ref_idx             Matrix of tx/rx pairing indices (in the order
%                       they're used in C).  If ignored, then all pairwise
%                       measurements are used (nTx x nRx).
%   angle_dims          Number of receiver angle of arrival dimensions
%                       reported (0 = no angles, 1 = azimuth, 2 = azimuth 
%                       & elevation)
%
% Outputs:
%   crlb    Lower bound on the error covariance matrix for an unbiased
%           PCL position estimator (Ndim x Ndim)
%
% Nicholas O'Donoughue
% 10 September 2021


if nargin < 9 || ~exist('angle_dims','var') || isempty(angle_dims)
    angle_dims = 0;
end
assert(angle_dims >=0 && angle_dims <= 2,'Error parsing angle dimensions command, must be between 0 and 2.');

do_doppler = ~isempty(v_tx) || ~isempty(v_rx) || ~isempty(v_src);

if nargin < 8 || ~exist('ref_idx','var')
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


%% Preprocess matrices and functions

% Construct Jacobian function handle
if do_doppler
    J = @(x,v) pcl.jacobian(x_tx,x_rx,x, v_tx, v_rx, v, ref_idx, angle_dims);
else
    J = @(x) pcl.jacobian(x_tx,x_rx,x,[],[],[],ref_idx, angle_dims);
end

% Preprocess covariance matrix
C_out = C;

% Pre-compute covariance matrix inverses
do_decomp = ~verLessThan('MATLAB','9.3');
if do_decomp
    % Starging in R2017b, MATLAB released the DECOMPOSITION function,
    % which can decompose matrices for faster computation of left- and
    % right-division in for loops.
    C_d = decomposition(C_out,'chol');
else
    % If DECOMPOSITION is unavailable, let's precompute the pseudo-inverse.
    C_inv = pinv(C_out);
end

% Initialize output variable
crlb = zeros([n_dim,n_dim,n_source]);

% Repeat CRLB for each of the n_source test positions
for idx =1:n_source
    this_x = x_src(:,idx);
    this_v = v_src(:,idx);

    % Evaluate Jacobian at x_i
    if do_doppler
        J_i = J(this_x, this_v);
    else
        J_i = J(this_x);
    end
    
    % Compute Fisher Information Matrix
    if do_decomp
        F = J_i*(C_d\J_i'); % Ndim x Ndim
    else
        F = J_i*C_inv*J_i';
    end
    
    if any(isnan(F(:))) || any(isinf(F(:)))
        % Problem is ill defined, Fisher Information Matrix cannot be
        % inverted
        crlb(:,:,idx) = NaN;
    else
        % Invert the Fisher Information Matrix to compute the CRLB
        crlb(:,:,idx) = pinv(F);
    end
end

