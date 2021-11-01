function J = grad_x(x_tdoa, x_source, tdoa_ref_idx, ~)
% J = grad_x(x_tdoa, x_source, tdoa_ref_idx, alpha_tdoa)          fdoa_ref_idx)
%
% Returns the gradient of hybrid measurements, with sensor uncertainties,
% with respect to target position, x, equation 6.23.
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

% Parse inputs
if nargin < 3 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

% Gradient with respect to target position is the same as the Jacobian,
% which was already derived in the absense of sensor position errors and
% measurement biases.  Just reference that function.
J = tdoa.jacobian(x_tdoa, x_source, tdoa_ref_idx);