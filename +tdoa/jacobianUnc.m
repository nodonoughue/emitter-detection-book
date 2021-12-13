function J = jacobianUnc(x_tdoa, x_source, tdoa_ref_idx, alpha_tdoa)
% J = jacobianUnc(x_tdoa, x_source, tdoa_ref_idx, alpha_tdoa);
%
% Returns the Jacobian matrix for hybrid set of AOA, TDOA, and FDOA
% measurements.
%
%
% INPUTS:
%   x_tdoa          TDOA sensor positions
%   x_source        Candidate source positions
%   tdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for TDOA measurements
%   alpha_tdoa
%
% OUTPUTS:
%   J               nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position
%
% Nicholas O'Donoughue
% 31 October 2021

% Parse Inputs
if nargin < 3 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

if nargin < 4 || ~exist('alpha_tdoa','var')
    alpha_tdoa = [];
end

% Gradient w.r.t source position
J_x = tdoa.grad_x(x_tdoa, x_source, tdoa_ref_idx, alpha_tdoa);

% Gradient w.r.t measurement biases
J_a = tdoa.grad_a(x_tdoa, x_source, tdoa_ref_idx, alpha_tdoa);

% Gradient w.r.t sensor position
J_b = tdoa.grad_b(x_tdoa, x_source, tdoa_ref_idx, alpha_tdoa);

% Combine component Jacobians
J = [J_x; J_a; J_b];