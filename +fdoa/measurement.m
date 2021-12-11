function rrdoa = measurement(x_sensor, v_sensor, x_source, ref_idx, alpha)
% rrdoa = measurement(x_sensor, v_sensor, x_source, ref_idx)
%
% Computed range rate difference measurements, using the
% final sensor as a common reference for all FDOA measurements.
%
% INPUTS:
%   x_sensor    nDim x nSensor array of sensor positions
%   v_sensor    nDim x nSensor array of sensor velocities
%   x_source    nDim x nSource array of source positions
%   ref_idx         Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings
%   alpha       nSensor x 1 array of range-rate bias terms
%
% OUTPUTS:
%   rrdoa       nSensor -1 x nSource array of RRDOA measurements [m/s]
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
[nDim1,nSensor1] = size(x_sensor);
[nDim2,nSensor2] = size(v_sensor);
[nDim3,nSource] = size(x_source);
if nDim1~=nDim2 || nSensor1 ~=nSensor2
    error('First two inputs must have macthing size');
end

if nDim1~=nDim3
    error('First dimension of all inputs must match');
end

if nargin < 3 || ~exist('ref_idx','var') || isempty(ref_idx)
    ref_idx = nSensor1;
end

% Parse Reference Sensor
n_sensor = size(x_sensor, 2);
[test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(ref_idx, n_sensor);

% Parse FDOA Bias
if nargin < 4 || ~exist('alpha','var') || isempty(alpha)
    rrdoa_bias = 0;
else
    rrdoa_bias = alpha(test_idx_vec) - alpha(ref_idx_vec);
end

% Compute distance from each source position to each sensor
dx = reshape(x_source,nDim1,1,nSource) - reshape(x_sensor,nDim1,nSensor1);
R = sqrt(sum(abs(dx).^2,1)); % 1 x nSensor x nSource

% Compute range rate from range and velocity
rr = reshape(sum(v_sensor.*dx./R,1),nSensor1,nSource); % nSensor x nSource

% Apply reference sensors to compute range rate difference for each sensor
% pair
rrdoa = rr(test_idx_vec,:) - rr(ref_idx_vec,:) + rrdoa_bias;