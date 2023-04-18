function crlb = computeCRLB(x_tx,x_rx,xs,C,ref_idx)
% crlb = computeCRLB(x_tdoa,xs,C,ref_idx)
%
% Computes the CRLB on position accuracy for source at location xs and
% sensors at locations in x1 (Ndim x N).  Ctdoa is an Nx1 vector of TOA
% variances at each of the N sensors.
%
% TODO: Incorporate range rate measurements
% TODO: Incorporate angle of arrival measurements
%
% Inputs:
%   
%   x_tx                Transmitter positions [m]
%   x_rx                Receiver positions [m]
%   xs                  (Ndim x M) array of source positions over which to 
%                       calculate CRLB
%   C                   Measurement Error Covariance Matrix [m^2]
%   ref_idx             Matrix of tx/rx pairing indices (in the order
%                       they're used in C).  If ignored, then all pairwise
%                       measurements are used (nTx x nRx)
%
% Outputs:
%   crlb    Lower bound on the error covariance matrix for an unbiased
%           PCL position estimator (Ndim x Ndim)
%
% Nicholas O'Donoughue
% 10 September 2021

% Parse inputs
if nargin < 5 || ~exist('ref_idx','var')
    ref_idx = [];
end
n_dim = size(x_tx,1);
n_source = size(xs,2);


% Construct Jacobian function handle
J = @(x) pcl.jacobian(x_tx,x_rx,x,[],[],[],ref_idx);

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
    this_x = xs(:,idx);
    
    % Evaluate Jacobian at x_i
    J_i = J(this_x);
    
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

