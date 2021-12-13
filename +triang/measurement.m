function psi = measurement(x_sensor, x_source, do2DAoA, alpha)
% psi = measurement(x_sensor, x_source, do2DAoA)
%
% Computes angle of arrival measurements.  If the input is 2D, then only
% the azimuth angle is report (atan(y/x)).  If the input is 3D, and the
% flag do2DAoA is missing or set to true, then the output contains two 
% sets of measurements, with azimuth reported first in psi(1:nSensor, :)
% and elevation reported second in psi(nSensor+1:end, :).
%
% INPUTS:
%   x_sensor    nDim x nSensor array of sensor positions
%   x_source    nDim x nSource array of source positions
%   do2DAoA     Boolean flag to activate 2D (az/el) AOA measurements
%               [default=True]
%   alpha       (Optional) nSensor x 1 vector of AOA biases (nSensor x 2 if
%               do2DAoA is true).  [default = 0]
%
% OUTPUTS:
%   psi         nSensor x nSource array of AOA measurements
%               (2*nSensor x nSource if 2D AOA measurements are generated)
%
% Nicholas O'Donoughue
% 1 July 2019

if nargin < 3 || isempty(do2DAoA)
    do2DAoA = true;
end

if nargin < 4 || isempty(alpha)
    alpha_az = 0;
    alpha_el = 0;
else
    alpha_az = alpha(:,1);
    if do2DAoA && size(alpha,2) > 1
        alpha_el = alpha(:,2);
    else
        alpha_el = 0;
    end
end

[nDim1,nSensor] = size(x_sensor);
[nDim2,nSource] = size(x_source);

if nDim1~=nDim2
    error('First dimension of all inputs must match');
end
if nDim1 < 2
    error('Must have at least two dimensions (x and y)');
end
nDim = nDim1;

dx = reshape(x_source,nDim,1,nSource) - reshape(x_sensor,nDim,nSensor);
        % nDim x nSensor x nSource
        
az = reshape(atan2(dx(2,:,:),dx(1,:,:)),nSensor,nSource) + alpha_az;

if nDim >2 && do2DAoA
    ground_rng = sqrt(sum(abs(dx(1:2,:,:)).^2,1));
    el = reshape(atan2(dx(3,:,:),ground_rng),nSensor,nSource) + alpha_el;
    
    psi = cat(1, az, el);
else
    psi = az;
end
