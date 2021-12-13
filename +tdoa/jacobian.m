function J = jacobian(x_sensor, x_source, ref_idx)
% J = jacobian(x_sensor, x_source, ref_idx)
%
% Returns the Jacobian matrix for TDOA of a source at x_source 
% (nDim x nSource) from sensors at x_sensor (nDim x nSensor).
%
% INPUTS:
%   x_sensor        nDim x nSensor vector of sensor positions
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

if nargin < 3 || ~exist('ref_idx','var') || isempty(ref_idx)
    ref_idx = nSensor;
end

% Parse Reference Sensors
n_sensor = size(x_sensor, 2);
[test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(ref_idx, n_sensor);

% Compute the range from each candidate source location to each sensor
dx = reshape(x_source,nDim,1,nSource) - reshape(x_sensor,nDim,nSensor);
R = sqrt(sum(abs(dx).^2,1)); % 1 x nSensor x nSource

% Remove any zero-range samples; replace with epsilon, to avoid a divide
% by zero error
R(R<1e-10) = 1e-10;

% Compute the Jacobians
J = dx(:,test_idx_vec,:)./R(:,test_idx_vec,:) ...
    - dx(:,ref_idx_vec,:)./R(:,ref_idx_vec,:);

