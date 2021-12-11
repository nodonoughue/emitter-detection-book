function rdoa = measurement(x_sensor, x_source, ref_idx, alpha)
% rdoa = measurement(x_sensor, x_source, ref_idx)
%
% Computes range difference measurements, using the
% final sensor as a common reference for all TDOA measurements.
%
% INPUTS:
%   x_sensor    nDim x nSensor array of sensor positions
%   x_source    nDim x nSource array of source positions
%   ref_idx     Either a scalar index for which sensor is the reference,
%               or a 2 x nPairing matrix of sensor pairing indices
%   alpha       nSensor x 1 vector of range bias terms
%
% OUTPUTS:
%   rdoa        nSensor -1 x nSource array of RDOA measurements
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
[nDim1,nSensor] = size(x_sensor);
[nDim2,nSource] = size(x_source);
if nDim1~=nDim2
    error('First dimension of all inputs must match');
end

if nargin < 3 || ~exist('ref_idx','var') || isempty(ref_idx)
    ref_idx = nSensor;
end

[test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(ref_idx, nSensor);

if nargin < 4 || ~exist('alpha','var') || isempty(alpha)
    rdoa_bias = 0;
else
    % Parse the TDOA bias
    rdoa_bias = (alpha(test_idx_vec) - alpha(ref_idx_vec));
end

% Compute range from each source to each sensor
dx = reshape(x_source,nDim1,1,nSource) - reshape(x_sensor,nDim1,nSensor);
R = reshape(sqrt(sum(abs(dx).^2,1)),nSensor,nSource); % nSensor x nSource

% Compute range difference for each pair of sensors
rdoa = R(test_idx_vec,:) - R(ref_idx_vec,:) + rdoa_bias;