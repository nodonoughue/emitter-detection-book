function [J,Jv] = jacobian(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source,tdoa_ref_idx,fdoa_ref_idx, v_source)
% [J,Jv] = jacobian(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source,tdoa_ref_idx,...
%                                                  fdoa_ref_idx, v_source)
%
% Returns the Jacobian matrix for hybrid set of AOA, TDOA, and FDOA
% measurements.
%
% If the target is moving, as specified by an optional fifth input
% v_source, then the Jacobian is provided with respect to both the target
% position and velocity.  This is only necessary if the geolocation
% algorithm is also solving for target velocity.  If target velocity is
% assumed known, or is not being estimated, then the source velocity can be
% subtracted from sensor velocity.
%
% INPUTS:
%   x_aoa           AOA sensor positions
%   x_tdoa          TDOA sensor positions
%   x_fdoa          FDOA sensor positions
%   v_fdoa          FDOA sensor velocities
%   x_source        Candidate source positions
%   tdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for TDOA measurements
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA measurements
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
% 1 July 2019

if nargin < 6 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

if nargin < 7 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

% Compute Jacobian for AOA measurements
if ~isempty(x_aoa)
    J_aoa = triang.jacobian(x_aoa, x_source);
else
    J_aoa = [];
end

% Compute Jacobian for TDOA measurements
if ~isempty(x_tdoa)
    J_tdoa= tdoa.jacobian(x_tdoa, x_source, tdoa_ref_idx);
else
    J_tdoa = [];
end

% Compute Jacobian for FDOA measurements
if ~isempty(x_fdoa) && ~isempty(v_fdoa)
    do_vel_jacobian = nargin > 4 && exist('v_source','var') && ~isempty(v_source);

    if do_vel_jacobian
        [J_fdoa, J_fdoa_v] = fdoa.jacobian(x_fdoa, v_fdoa, x_source, fdoa_ref_idx, v_source);
    else
        J_fdoa = fdoa.jacobian(x_fdoa, v_fdoa, x_source, fdoa_ref_idx);
        J_fdoa_v = [];
    end
else
    do_vel_jacobian = false;

    J_fdoa = [];
    J_fdoa_v = [];
end

%% Combine component Jacobians
J = [J_aoa, J_tdoa, J_fdoa];

if do_vel_jacobian
    num_rows = size(J_fdoa_v,1);
    num_zero_cols = size(J_aoa,2) + size(J_tdoa,2);
    Jv = cat(2,zeros(num_rows, num_zero_cols), J_fdoa_v);
else
    Jv = [];
end