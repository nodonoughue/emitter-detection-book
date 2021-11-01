function J = grad_b(x_tdoa, x_source, tdoa_ref_idx, ~)
% J = grad_b(x_tdoa, x_source,tdoa_ref_idx,alpha_tdoa)
%
% Returns the gradient of hybrid measurements, with sensor uncertainties,
% with respect to sensor position and velocity. Equation 6.27
%
%
% INPUTS:
%   x_tdoa          TDOA sensor positions
%   x_source        Candidate source positions
%   tdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for TDOA measurements
%   alpha_tdoa
%
% OUTPUTS:
%   J               nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position
%
% Nicholas O'Donoughue
% 31 October 2021

n_dim = size(x_tdoa,1);
n_tdoa = size(x_tdoa,2);
n_source = size(x_source,2);

if nargin < 3 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

% Compute Pointing Vectors
dx = x_tdoa - reshape(x_source,n_dim,1,n_source); % nDim x nSensor x nSource
Rn = sqrt(sum(abs(dx).^2,1)); % Euclidean norm for each offset vector
dx_norm = dx ./ Rn;

% Parse reference index vector
[test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(tdoa_ref_idx, n_tdoa);

% Build Gradient
n_msmt = numel(test_idx_vec);
J = zeros(n_dim*n_tdoa, n_msmt, n_source);
for idx=1:n_msmt
    this_test = test_idx_vec(idx);
    this_ref = ref_idx_vec(idx);
    
    J((1:n_dim) + n_dim*(this_test-1),idx,:) = -dx_norm(:,this_test,:);
    J((1:n_dim) + n_dim*(this_ref-1),idx,:)  =  dx_norm(:,this_ref,:);
end
