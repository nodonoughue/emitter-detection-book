function J = grad_b(x_aoa, x_source, do2DAoA, alpha_aoa)
% J = grad_b(x_aoa, x_source, do2DAoA, alpha_aoa)
%
% Returns the gradient of hybrid measurements, with sensor uncertainties,
% with respect to sensor position and velocity. Equation 6.17.
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

if nargin < 3 || ~exist('do2DAoA','var')
    do2DAoA = size(x_source,1)==3;
end

if nargin < 4 || ~exist('alpha_aoa','var')
    alpha_aoa = [];
end

n_dim = size(x_aoa, 1);
n_aoa = size(x_aoa, 2);
n_source = size(x_source,2);

% Compute Jacobian for Azimuth measurements
% Equation 6.18 and 6.19 show that the gradient
% w.r.t. sensor position is the negative of the gradient w.r.t. target
% position, but resampled to be block diagonal.
J_x = triang.grad_x(x_aoa, x_source, do2DAoA, alpha_aoa);
        % Result is 3 x n_aoa or 3 x (2*n_aoa)

J = zeros(n_dim*n_aoa, n_aoa, n_source);
for idx=1:n_aoa
    J((1:n_dim)+n_dim*(idx-1),idx,:) = J_x(:,idx,:);
end

% Compute Jacobian for Elevation measurements
if do2DAoA && n_dim == 3
    J_b_el = zeros(3*n_aoa, n_aoa,n_source);
    for idx=1:n_aoa
        J_b_el((1:n_dim)+n_dim*(idx-1),idx,:) = J_x(:,idx+n_aoa,:);
    end
    
    J = cat(2,J, J_b_el);
end
