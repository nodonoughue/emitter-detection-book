function J = jacobian(x_sensor,v_sensor,x_source, ref_idx)
% J = jacobian(x_sensor,v_sensor,x_source, ref_idx)
%
% Returns the Jacobian matrix for FDOA of a source at x_source 
% (nDim x nSource) from sensors at x_sensor (nDim x nSensor) with velocity 
% v_sensor.
%
% INPUTS:
%   x_sensor        nDim x nSensor vector of sensor positions
%   v_sensor        nDim x nSensor vector of sensor velocities
%   x_source        nDim x nSource vector of source positions
%   ref_idx         Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings
%
% OUTPUTS:
%   J               nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
[nDim1,nSensor] = size(x_sensor);
[nDim2,nSource] = size(x_source);

if nDim1 ~= nDim2
    error('Input variables must match along first dimension.');
end
nDim = nDim1;

if nargin < 4 || ~exist('ref_idx','var') || isempty(ref_idx)
    ref_idx = nSensor;
end

if isscalar(ref_idx)
    test_idx_vec = setdiff(1:nSensor,ref_idx);
    ref_idx_vec = ref_idx;
else
    test_idx_vec = ref_idx(1,:);
    ref_idx_vec = ref_idx(2,:);
end

%% Compute the Offset Vectors
dx = x_sensor - reshape(x_source,nDim,1,nSource); % nDim x nSensor x nSource
Rn = sqrt(sum(abs(dx).^2,1)); % Euclidean norm for each offset vector
dx_norm = dx ./ sqrt(sum(abs(dx).^2,1));
Px = reshape(dx_norm,nDim,1,nSensor,nSource).*reshape(dx_norm,1,nDim,nSensor,nSource); % nDim x nDim x nSensor

%% Compute the gradient of R_n
nabla_Rn = squeeze(sum((eye(nDim) - Px) ...
                .* reshape(v_sensor./Rn,1,nDim,nSensor),2));
      % nDim x nSensor x nSource

%% Take the reference of each w.r.t. to the N-th
J = nabla_Rn(:,test_idx_vec,:) - nabla_Rn(:,ref_idx_vec,:);
    % nDim x nPair x nSource
    
    
