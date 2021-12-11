function J = jacobianUnc(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source,tdoa_ref_idx,fdoa_ref_idx, do2DAoA, alpha_aoa, alpha_tdoa, alpha_fdoa)
% J = jacobianUnc(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source,tdoa_ref_idx,
%                 fdoa_ref_idx, do2DAoA, alpha_aoa, alpha_tdoa, alpha_fdoa)
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
%   do2DAOA
%   alpha_aoa
%   alpha_tdoa
%   alpha_fdoa
%
% OUTPUTS:
%   J               nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position
%
% Nicholas O'Donoughue
% 31 October 2021

% Parse Inputs
if nargin < 6 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

if nargin < 7 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

if nargin < 8 || ~exist('do2DAoA','var')
    do2DAoA = [];
end

if nargin < 9 || ~exist('alpha_aoa','var')
    alpha_aoa = [];
end

if nargin < 10 || ~exist('alpha_tdoa','var')
    alpha_tdoa = [];
end

if nargin < 11 || ~exist('alpha_fdoa','var')
    alpha_fdoa = [];
end

% Gradient w.r.t source position
J_x = hybrid.grad_x(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source,tdoa_ref_idx,fdoa_ref_idx, do2DAoA, alpha_aoa, alpha_tdoa, alpha_fdoa);

% Gradient w.r.t measurement biases
J_a = hybrid.grad_a(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source,tdoa_ref_idx,fdoa_ref_idx, do2DAoA, alpha_aoa, alpha_tdoa, alpha_fdoa);

% Gradient w.r.t sensor position
J_b = hybrid.grad_b(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source,tdoa_ref_idx,fdoa_ref_idx, do2DAoA, alpha_aoa, alpha_tdoa, alpha_fdoa);

% Combine component Jacobians
J = [J_x; J_a; J_b];