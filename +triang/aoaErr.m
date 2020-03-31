function [epsilon,x_vec,y_vec] = aoaErr(x_sensor,C,x_source,x_max,numPts)
% [epsilon,x_vec,y_vec] = aoaErr(x_sensor,C,x0,x_max,numPts)
%
% Construct a 2-D field from -x_max to +x_max, using numPts in each
% dimension.  For each point, compute the AOA solution for each sensor, and 
% compare to the AOA solution from the true emitter position.
%
% INPUTS:
%   x_sensor    nDim x N matrix of sensor positions
%   x_source    nDim x 1 matrix of true emitter position
%   x_max       nDim x 1 (or scalar) vector of maximum offset from origin
%               for plotting
%   numPts      Number of test points along each dimension
%
% OUTPUTS:
%   epsilon     2-D plot of AOA error
%   x_vec       Vector of x coordinates for 2-D plot
%   y_vec       Vector of y coordinates for 2-D plot
%
% Nicholas O'Donoughue
% 20 March 2020

%% Compute the True AOA measurement
% Compute true range difference measurements; default condition is to use
% the final sensor as the reference for all difference measurements.
psi = triang.measurement(x_sensor,x_source);

% Preprocess covariance matrix inverses
C_out = C(1:end-1,1:end-1) + C(end,end);
C_d = decomposition(C_out);
    
%% Set up test points
xx_vec = x_max(:).*repmat(linspace(-1,1,numPts),2,1);
x_vec = xx_vec(1,:);
y_vec = xx_vec(2,:);
[XX,YY] = ndgrid(xx_vec(1,:),xx_vec(2,:));
x_plot = [XX(:),YY(:)]'; % 2 x numPts^2

epsilon = zeros(size(XX));
for idx_pt = 1:numel(XX)
    x_i = x_plot(:,idx_pt);
    psi_i = triang.measurement(x_sensor,x_i);
    err = psi-psi_i;
    epsilon(idx_pt) = (err'/C_d*err);
end
