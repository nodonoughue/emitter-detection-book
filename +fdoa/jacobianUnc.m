function J = jacobianUnc(x_fdoa, v_fdoa, x_source, fdoa_ref_idx, alpha_fdoa)
% J = jacobianUnc(x_fdoa, v_fdoa, x_source, fdoa_ref_idx,  alpha_fdoa)
%
% Returns the Jacobian matrix for hybrid set of AOA, TDOA, and FDOA
% measurements.
%
%
% INPUTS:
%   x_fdoa          FDOA sensor positions
%   v_fdoa          FDOA sensor velocities
%   x_source        Candidate source positions
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA measurements
%   alpha_fdoa
%
% OUTPUTS:
%   J               nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse Inputs
if nargin < 4 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

if nargin < 5 || ~exist('alpha_fdoa','var')
    alpha_fdoa = [];
end

% Gradient w.r.t source position
J_x = fdoa.grad_x(x_fdoa, v_fdoa, x_source, fdoa_ref_idx, alpha_fdoa);

% Gradient w.r.t measurement biases
J_a = fdoa.grad_a(x_fdoa, v_fdoa, x_source, fdoa_ref_idx, alpha_fdoa);

% Gradient w.r.t sensor position
J_b = fdoa.grad_b(x_fdoa, v_fdoa, x_source, fdoa_ref_idx, alpha_fdoa);

% Combine component Jacobians
J = [J_x; J_a; J_b];