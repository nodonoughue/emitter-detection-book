function J = jacobian(x_rdr, x_tgt, v_rdr, v_tgt, angle_dims)
% J = jacobian(x_tx, x_tgt, v_tx, v_tgt, angle_dims)
%
% Returns the Jacobian matrix for monostatic range, range-rate, and angle
% measurements of a source at x_source (nDim x nSource).
%
% INPUTS:
%   x_rdr               Radar positions [m]
%   x_tgt               Target positions [m]
%   v_rdr               Radar velocities [m/s]
%   v_tgt               Target velocities [m/s]
%   angle_dims       Flag determining how many angle of arrival dimensions 
%                    to report (0 =  no angle data, 1 = azimuth, 2 =
%                    azimuth & elevation)
%
% OUTPUTS:
%   J               nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position
%
% Nicholas O'Donoughue
% 15 September 2022

%% Parse inputs
[nDim,nRdr] = size(x_rdr);
[nDim3,nTgt] = size(x_tgt);

assert(nDim == nDim2 && nDim2 == nDim3, ...
       'Input variables must match along first dimension.');

if nargin < 5 || ~exist('angle_dims','var') || isempty(angle_dims)
    angle_dims = 0;
end
assert(angle_dims >=0 && angle_dims <= 2,'Error parsing angle dimensions command, must be between 0 and 2.');

do_doppler = nargin >= 3 && (~isempty(v_rdr) || ~isempty(v_tgt));
if do_doppler && isempty(v_rdr)
    v_rdr = zeros(size(x_rdr));
end

if do_doppler && isempty(v_tgt)
    v_tgt = zeros(size(x_tgt));
end

%% Range Jacobian
epsilon = 1e-15; % minimum range; to avoid divide by zero errors

% Compute Tx-Tgt Ranges
dx = reshape(x_tgt,nDim,nTgt) - reshape(x_rdr,nDim,nRdr);
range = max(sqrt(sum(abs(dx).^2,1)),epsilon); % 1 x nTx x nTgt
J_rng = dx./range; % Tx-tgt component of jacobian

%% Range-Rate Jacobian
if do_doppler
    u_tx = J_rng; % The transmitter-tgt jacobian is u_m(x) for transmitter m
    
    % Projection matrix
    Proj_ortho = eye(nDim) - reshape(u_tx,nDim,1,nRdr,nTgt).*reshape(u_tx,1,nDim,nRdr,nTgt);
    
    % Scaled velocity vector
    uv = (reshape(v_rdr,nDim,nRdr) - reshape(v_tgt,nDim,1,nTgt)) ./ range; % nDim x nRdr x nTgt
    J_velrng = squeeze(sum(Proj_ortho.*reshape(uv,1,nDim,nRdr,nTgt),2)); % nDim x nRdr x nTgt

    % Jacobian of range-rate w.r.t velocity is equivalent to the negative 
    % of the jacobian of range w.r.t. position
    J_vel = -J_rng;
end

%% Angle of Arrival Measurements
if angle_dims > 0 
    do2DAOA = (angle_dims==2);
    
    J_psi = triang.jacobian(x_rdr, x_tgt, do2DAOA);
else
    J_psi = [];
end

J_rng_full = cat(2,J_rng,J_velrng,J_psi); % Jacobian w.r.t. position
J_vel_full = cat(2,zeros(size(J_rng)), J_vel, zeros(size(J_psi))); % Jacobian w.r.t. velocity
J = cat(1,J_rng_full,J_vel_full); % Stack them for the full 6-state Jacobian