function crlb = computeCRLB(x_fdoa,v_fdoa,xs,C,ref_idx)
% crlb = computeCRLB(x_fdoa,v_fdoa,xs,C,ref_idx)
%
% Computes the CRLB on position accuracy for source at location xs and
% sensors at locations in x_fdoa (Ndim x N) with velocity v_fdoa.  
% C is an Nx1 vector of FOA variances at each of the N sensors, and ref_idx
% defines the reference sensor(s) used for FDOA.
%
% Inputs:
%   x_fdoa      (Ndim x N) array of FDOA sensor positions
%   v_fdoa      (Ndim x N) array of FDOA sensor velocities
%   xs          (Ndim x M) array of source positions over which to 
%               calculate CRLB
%   C           Covariance matrix for frequency estimates at the N
%               FDOA sensors [Hz^2]
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
n_dim = size(x_fdoa,1);
n_source = size(xs,2);

% Define the Jacobian function
J = @(x) fdoa.jacobian(x_fdoa,v_fdoa,x,ref_idx);

% Pre-process the covariance matrix for faster computation
C_d = decomposition(C);

% Initialize output variable
crlb = zeros([n_dim,n_dim,n_source]);

% Repeat CRLB for each of the n_source test positions
for idx =1:n_source
    this_x = xs(:,idx);
    J_i = J(this_x);
    F = J_i/C_d*J_i'; % Ndim x Ndim
    if any(isnan(F(:))) || any(isinf(F(:)))
        % Problem is ill defined, Fisher Information Matrix cannot be
        % inverted
        crlb(:,:,idx) = NaN;
        continue;
    end
    
    finv = pinv(F);
        
    % Condition the errors; set negative values to 0 for now
    crlb(:,:,idx) = finv;
end