function crlb = computeCRLB(x_rdr,x_tgt,v_rdr,v_tgt,C,angle_dims,print_status)
% crlb = computeCRLB(x_rdr,x_tgt,v_rdr,v_tgt,C,angle_dims)
%
% Computes the CRLB on position accuracy for target at location x_tgt and
% velocity v_tgt, given radars at locations x_rdr with velocities v_rdr.
% Assumes perfect knowledge of sensor and transmitter parameters.
%
% Inputs:
%   
%   x_rdr               Radar positions [m]
%   x_tgt               Candidate target position [m]
%   v_rdr               Radar velocities [m/s]
%   v_tgt               Candidate target velocity [m/s]
%   C                   Measurement Error Covariance Matrix [m^2, m^2/s^2,
%                       rad^2]
%   angle_dims          Number of receiver angle of arrival dimensions
%                       reported (0 = no angles, 1 = azimuth, 2 = azimuth 
%                       & elevation)
%   print_status        Optional flag, controlling print statement on
%                       how many solutions were nearly singular.
%                       Default=True
%
% Outputs:
%   crlb    Lower bound on the error covariance matrix for an unbiased
%           PCL position estimator (Ndim x Ndim)
%
% Nicholas O'Donoughue
% 15 September 2022

%% Turn off warnings for badly scaled matrices
warning('off','MATLAB:nearlySingularMatrix');
singular_matrix_errors = 0;

%% Parse inputs
if nargin < 7 || ~exist('print_status','var') || isempty(print_status)
    print_status = 0;
end

if nargin < 6 || ~exist('angle_dims','var') || isempty(angle_dims)
    angle_dims = 0;
end
assert(angle_dims >=0 && angle_dims <= 2,'Error parsing angle dimensions command, must be between 0 and 2.');

do_doppler = ~isempty(v_rdr) || ~isempty(v_tgt);

if do_doppler && isempty(v_rdr)
    v_rdr = zeros(size(x_rdr));
end

if do_doppler && isempty(v_tgt)
    v_tgt = zeros(size(x_tgt));
end

% Handle position-specific CRLB inputs
if numel(size(C)) > 2
    assert(size(C,3)==size(x_tgt,2),'For target position specific covariance matrices, the third dimension of the covariance matrix must match the number of target positions.');
    % Make sure the velocity is similarly sized
    if size(v_tgt,2)==1 || size(v_tgt,1)==1
        v_tgt = v_tgt(:) * ones(1,size(x_tgt,2));
    end
    if do_doppler
        crlb_cell = arrayfun(@(idx) pcl.computeCRLB(x_rdr,x_tgt(:,idx),v_rdr,v_tgt(:,idx),C(:,:,idx),angle_dims,false),1:size(C,3),'UniformOutput',false);
    else
        crlb_cell = arrayfun(@(idx) pcl.computeCRLB(x_rdr,x_tgt(:,idx),[],[],C(:,:,idx),angle_dims,false),1:size(C,3),'UniformOutput',false);
    end
        % n_tgt x 1 cell array of CRLBs, reshape
    crlb = cell2mat(reshape(crlb_cell,1,1,size(C,3)));
    return;
end

% Check source pos/vel inputs
if do_doppler
    size_pos = size(x_tgt);
    size_vel = size(v_tgt);
    assert( ~any((size_pos ~= size_vel) & (size_pos > 1) & (size_vel > 1)), 'Non-singleton dimensions of source position and velocity must match.');

    in_size = max(size_pos,size_vel);
    out_size = in_size(2:end); % drop first dimension
    n_dim = in_size(1);
    n_source_pos = prod(out_size);

    % Broadcast x_source and v_source to common size
    x_tgt = reshape(x_tgt.*ones(in_size), n_dim, n_source_pos);
    v_tgt = reshape(v_tgt.*ones(in_size), n_dim, n_source_pos);
else
    in_size = size(x_tgt);
    out_size = in_size(2:end); % drop first dimension
    n_dim = in_size(1);
    n_source_pos = prod(out_size);

    x_tgt = reshape(x_tgt,n_dim,n_source_pos);
    v_tgt = zeros(n_dim,n_source_pos);
end


%% Preprocess matrices and functions

% Construct Jacobian function handle
if do_doppler
    J = @(x,v) radar.jacobian(x_rdr,x, v_rdr, v, angle_dims);
else
    J = @(x) radar.jacobian(x_rdr,x,[],[], angle_dims);
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
    % If DECOMPOSITION is unavailable, let's precompute the inverse.
    C_inv = inv(C_out);
end

% Initialize output variable
if do_doppler
    crlb = zeros([2*n_dim,2*n_dim,n_source_pos]);
else
    crlb = zeros([n_dim,n_dim,n_source_pos]);
end
% Repeat CRLB for each of the n_source test positions
for idx =1:n_source_pos
    this_x = x_tgt(:,idx);
    this_v = v_tgt(:,idx);

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
        crlb(:,:,idx) = inv(F);
        
        % Check for singular matrices and count
        ME = warning('query','last');
        if ~isempty(ME) && strcmpi(ME.identifier,'MATLAB:nearlySingularMatrix')
            singular_matrix_errors = singular_matrix_errors + 1;
        end
        % Clear the last warning, so we don't double count it
        lastwarn('');
    end
end

warning('on','MATLAB:nearlySingularMatrix');
if singular_matrix_errors > 0 && print_status
    fprintf('computeCRLB encountered %d degenerate cases out of %d source positions. Problem potentially ill-defined.\n',singular_matrix_errors, n_source_pos);
end