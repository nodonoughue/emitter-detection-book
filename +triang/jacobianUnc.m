function J = jacobianUnc(x_aoa, x_source, do2DAoA, alpha_aoa)
% J = jacobianUnc(x_aoa, x_source, do2DAoA, alpha_aoa)
%
% Returns the Jacobian matrix for hybrid set of AOA, TDOA, and FDOA
% measurements.
%
%
% INPUTS:
%   x_aoa           AOA sensor positions
%   x_source        Candidate source positions
%   do2DAOA
%   alpha_aoa
%
% OUTPUTS:
%   J               nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position
%
% Nicholas O'Donoughue
% 31 October 2021

% Parse Inputs
if nargin < 3 || ~exist('do2DAoA','var')
    do2DAoA = [];
end

if nargin < 4 || ~exist('alpha_aoa','var')
    alpha_aoa = [];
end

% Gradient w.r.t source position
J_x = triang.grad_x(x_aoa, x_source, do2DAoA, alpha_aoa);

% Gradient w.r.t measurement biases
J_a = triang.grad_a(x_aoa, x_source, do2DAoA, alpha_aoa);

% Gradient w.r.t sensor position
J_b = triang.grad_b(x_aoa, x_source, do2DAoA, alpha_aoa);

% Combine component Jacobians
J = [J_x; J_a; J_b];