function J = grad_x(x_fdoa, v_fdoa, x_source, fdoa_ref_idx, ~)
% J = grad_x(x_fdoa, v_fdoa, x_source, fdoa_ref_idx, alpha_fdoa)
%
% Returns the gradient of hybrid measurements, with sensor uncertainties,
% with respect to target position, x. Equation 6.31.
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
% 31 October 2021

% Parse Inputs
if nargin < 4 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end


% The gradient with respect to sensor position x is the same as the
% Jacobian that was previously implemented without sensor position errors
% or measurement biases, so just reference that function.
J = fdoa.jacobian(x_fdoa, v_fdoa, x_source, fdoa_ref_idx);
