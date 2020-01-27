function J = jacobian(x_sensor, x_source)
% J = jacobian(x_sensor, x_source)
%
% Returns the Jacobian matrix for triangulation of a source at x_source 
% (nDim x nSource) from sensors at x_sensor (nDim x nSensor)
%
% INPUTS:
%   x_sensor        nDim x nSensor vector of sensor positions
%   x_source        nDim x nSource vector of source positions
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
if nDim1~=2
    error('Jacobian not currently defined for 3-D angle of arrival.');
end

% Compute offset/range from each sensor to each source position
dx = reshape(x_source,nDim1,1,nSource) - reshape(x_sensor,nDim1,nSensor);
R2 = sum(abs(dx).^2,1); % 1 x nSensor x nSource

% Build Jacobian
J = cat(1,-dx(2,:,:),dx(1,:,:))./R2;