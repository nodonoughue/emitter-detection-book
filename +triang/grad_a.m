function J = grad_a(x_aoa, x_source, do2DAoA, ~)
% J = grad_a(x_aoa, x_source, do2DAoA, alpha_aoa)
%
% Returns the gradient of hybrid measurements, with sensor uncertainties,
% with respect to sensor measurement bias, alpha. Equation 6.20
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

% Compute Jacobian for AOA measurements
n_dim = size(x_aoa, 1);
n_aoa = size(x_aoa, 2);
n_source = size(x_source,2);
if do2DAoA && n_dim ==3
    n_out = 2*n_aoa;
else
    n_out = n_aoa;
end

J = repmat(eye(n_out),1,1,n_source);