function [psi, psi_elev] = measurement(x_sensor, x_source)
% psi, psi_elev = measurement(x_sensor, x_source)
%
% Computes angle of arrival measurements.
%
% INPUTS:
%   x_sensor    nDim x nSensor array of sensor positions
%   x_source    nDim x nSource array of source positions
%
% OUTPUTS:
%   psi         nSensor x nSource array of azimuth measurements (2D or 3D)
%   psi_elev    nSensor x nSource array of elevation measurements (3D)
%
% Nicholas O'Donoughue
% 1 July 2019

[nDim1,nSensor] = size(x_sensor);
[nDim2,nSource] = size(x_source);

if nDim1~=nDim2
    error('First dimension of all inputs must match');
end
nDim = nDim1;

dx = reshape(x_source,nDim,1,nSource) - reshape(x_sensor,nDim,nSensor);
        % nDim x nSensor x nSource
        
psi = reshape(atan2(dx(2,:,:),dx(1,:,:)),nSensor,nSource);

if nDim1 == 3
    % Compute elevation angle
    psi_elev = reshape(atan2(dx(3,:,:),sqrt(dx(1,:,:).^2+dx(2,:,:).^2)),...
                        nSensor,nSource);
else
    psi_elev = [];
end