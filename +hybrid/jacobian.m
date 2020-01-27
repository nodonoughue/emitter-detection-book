function J = jacobian(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source,tdoa_ref_idx,fdoa_ref_idx)
% J = jacobian(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source,tdoa_ref_idx,...
%                                                      fdoa_ref_idx)
%
% Returns the Jacobian matrix for hybrid set of AOA, TDOA, and FDOA
% measurements.
%
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
%
% OUTPUTS:
%   J               nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position
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
    J_fdoa= fdoa.jacobian(x_fdoa, v_fdoa, x_source, fdoa_ref_idx);
else
    J_fdoa = [];
end

% Combine component Jacobians
J = [J_aoa, J_tdoa, J_fdoa];