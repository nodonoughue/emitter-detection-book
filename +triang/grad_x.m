function J = grad_x(x_aoa, x_source, do2DAoA, ~)
% J = grad_x(x_aoa, x_source, do2DAoA, alpha_aoa)
%
% Returns the gradient of hybrid measurements, with sensor uncertainties,
% with respect to target position, x. Equation 6.14.
%
%
% INPUTS:
%   x_aoa           AOA sensor positions
%   x_source        Candidate source positions
%   do2DAoA
%   alpha_aoa
%
% OUTPUTS:
%   J               nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position
%
% Nicholas O'Donoughue
% 31 October 2021


% Parse inputs
if nargin < 3 || ~exist('do2DAoA','var')
    do2DAoA = size(x_source,1)==3;
end


% The gradient is unchanged; which we have already implemented in the
% absence of sensor position error and msmt bias, in triang.jacobian.
% So, we'll just reference that.
J = triang.jacobian(x_aoa, x_source, do2DAoA);