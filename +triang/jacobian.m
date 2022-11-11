function J = jacobian(x_sensor, x_source, do2DAoA)
% J = jacobian(x_sensor, x_source, do2DAoA)
%
% Returns the Jacobian matrix for triangulation of a source at x_source 
% (nDim x nSource) from sensors at x_sensor (nDim x nSensor)
%
% INPUTS:
%   x_sensor        nDim x nSensor vector of sensor positions
%   x_source        nDim x nSource vector of source positions
%   do2DAoA         Boolean flag indicating whether 2D AOA (azimuth and
%                   elevation) is to be used
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

if nargin < 3 || ~exist('do2DAoA','var')
    do2DAoA = true;
end

% Compute offset/range from each sensor to each source position
dx = reshape(x_source,nDim1,1,nSource) - reshape(x_sensor,nDim1,nSensor);
R2 = sum(abs(dx).^2,1); % 1 x nSensor x nSource  -- slant range
GR2 = sum(abs(dx(1:2,:,:)).^2,1); % 1 x nSensor x nSource -- ground range

% Build Jacobian
J = cat(1,-dx(2,:,:),dx(1,:,:))./GR2;

% Elevation Angle Jacobian
if do2DAoA && nDim1 == 3

    % Add a z dimension to the azimuth Jacobian
    J = cat(1, J, zeros(1, nSensor, nSource));
    
    % Jacobian for elevation
    % del_x(phi_n(x)) = -(x-xn)(z-zn)/ ground_range * slant_range^2
    % del_y(phi_n(x)) = -(y-yn)(z-zn)/ ground_range * slant_range^2
    % del_z(phi_n(x)) = ground_range / slant_range^2
    J_el = cat(1,-dx(1,:,:).*dx(3,:,:)./sqrt(GR2), ...
                 -dx(2,:,:).*dx(3,:,:)./sqrt(GR2), ...
                 sqrt(GR2))./R2;
             
    
    % The elevation measurements are concatenated after the azimuth
    % measurements, so the Jacobian is concatenated in the second (nSensor)
    % dimension
    J = cat(2, J, J_el);
elseif nDim1 == 3
    % It's a 3D/az-only problem
    J = cat(1,J, zeros(1,size(J,2)));
end
