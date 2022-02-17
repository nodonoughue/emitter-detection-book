function [J, Jv] = jacobian(x_sensor,v_sensor,x_source, ref_idx, v_source)
% [J, Jv] = jacobian(x_sensor,v_sensor,x_source, ref_idx, v_source)
%
% Returns the Jacobian matrix for FDOA of a source at x_source 
% (nDim x nSource) from sensors at x_sensor (nDim x nSensor) with velocity 
% v_sensor.
%
% If the target is moving, as specified by an optional fifth input
% v_source, then the Jacobian is provided with respect to both the target
% position and velocity.  This is only necessary if the geolocation
% algorithm is also solving for target velocity.  If target velocity is
% assumed known, or is not being estimated, then the source velocity can be
% subtracted from sensor velocity.
%
% INPUTS:
%   x_sensor        nDim x nSensor vector of sensor positions
%   v_sensor        nDim x nSensor vector of sensor velocities
%   x_source        nDim x nSource vector of source positions
%   ref_idx         Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings
%   v_source        [Optional] nDim x nSource vector of source velocities
%                   Target assumed stationary if not provided.
%
% OUTPUTS:
%   J               nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position
%   Jv              nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position (if v_source is
%                   provided).
%
% Nicholas O'Donoughue
% 24 January 2022

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

do_vel_jacobian = nargin > 4 && exist('v_source','var') && ~isempty(v_source);    

% Parse Reference Sensor
n_sensor = size(x_sensor, 2);
[test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(ref_idx, n_sensor);

%% Compute the Offset Vectors
dx = x_sensor - reshape(x_source,nDim,1,nSource); % nDim x nSensor x nSource
Rn = sqrt(sum(abs(dx).^2,1)); % Euclidean norm for each offset vector
dx_norm = dx ./ sqrt(sum(abs(dx).^2,1));
Px = reshape(dx_norm,nDim,1,nSensor,nSource).*reshape(dx_norm,1,nDim,nSensor,nSource); % nDim x nDim x nSensor

if do_vel_jacobian
    % There's a non-zero velocity
    dv = v_sensor - v_source;
else
    dv = v_sensor;
end

%% Compute the gradient of R_n
nabla_Rn = squeeze(sum((eye(nDim) - Px) ...
                .* reshape(dv./Rn,1,nDim,nSensor),2));
      % nDim x nSensor x nSource

%% Take the reference of each w.r.t. to the N-th
J = nabla_Rn(:,test_idx_vec,:) - nabla_Rn(:,ref_idx_vec,:);
    % nDim x nPair x nSource
    
%% Velocity Jacobian
if do_vel_jacobian
    Jv = dx_norm(:,test_idx_vec,:) - dx_norm(:,ref_idx_vec,:);
else
    Jv = [];
end
