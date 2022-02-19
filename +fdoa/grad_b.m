function J = grad_b(x_fdoa, v_fdoa, x_source, fdoa_ref_idx, ~)
% J = grad_b(x_fdoa, v_fdoa, x_source, fdoa_ref_idx, alpha_fdoa)
%
% Returns the gradient of hybrid measurements, with sensor uncertainties,
% with respect to sensor position and velocity.
%
%
% INPUTS:
%   x_fdoa          FDOA sensor positions
%   v_fdoa          FDOA sensor velocities
%   x_source        Candidate source positions
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA measurements
%   alpha_fdoa
%
% OUTPUTS:
%   J               nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position
%
% Nicholas O'Donoughue
% 31 October 2021

if nargin < 4 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

n_dim = size(x_fdoa,1);
n_fdoa = size(x_fdoa,2);
n_source = size(x_source,2);

if nargin < 4 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

% Compute Pointing Vectors
dx = x_fdoa - reshape(x_source,n_dim,1,n_source); % nDim x nSensor x nSource
Rn = sqrt(sum(abs(dx).^2,1)); % Euclidean norm for each offset vector
dx_norm = dx ./ Rn;

Px = reshape(dx_norm,n_dim,1,n_fdoa,n_source).*reshape(conj(dx_norm),1,n_dim,n_fdoa,n_source);

%% Compute the gradient of R_n
nabla_Rn = squeeze(sum((eye(n_dim) - Px) ...
                .* reshape(v_fdoa./Rn,1,n_dim,n_fdoa,n_source),2));
      % nDim x nSensor x nSource

%% Parse reference index vector
[test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(fdoa_ref_idx, n_fdoa);

% Build Gradient
n_msmt = numel(test_idx_vec);
Jx = zeros(n_dim*n_fdoa, n_msmt, n_source);
Jv = zeros(n_dim*n_fdoa, n_msmt, n_source);
for idx=1:n_msmt
    this_test = test_idx_vec(idx);
    this_ref = ref_idx_vec(idx);
    
    % Gradient w.r.t. Sensor Position (eq 6.35)
    Jx((1:n_dim) + n_dim*(this_test-1),idx,:) =  nabla_Rn(:,this_test,:);
    Jx((1:n_dim) + n_dim*(this_ref-1),idx,:)  = -nabla_Rn(:,this_ref,:);

    % Gradient w.r.t Sensor Velocity (eq 6.36)
    Jv((1:n_dim) + n_dim*(this_test-1),idx,:) = -dx_norm(:,this_test,:);
    Jv((1:n_dim) + n_dim*(this_ref-1),idx,:)  =  dx_norm(:,this_ref,:);
    
end

%% Combine the Jacobian
% eq 6.34
J = cat(1,Jx,Jv);
