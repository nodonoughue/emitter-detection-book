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

% Make sure the covariance matrix is invertible
C_out = utils.ensureInvertible(C);

do_decomp = ~verLessThan('MATLAB','9.3');
if do_decomp
    % Starging in R2017b, MATLAB released the DECOMPOSITION function,
    % which can decompose matrices for faster computation of left- and
    % right-division in for loops.
    C_d = decomposition(C_out,'chol');
else
    % If DECOMPOSITION is unavailable, let's precompute the pseudo-inverse.
    C_inv = pinv(C_out);
end
    
%% Set up test points
xx_vec = x_max(:).*repmat(linspace(-1,1,numPts),2,1);
x_vec = xx_vec(1,:);
y_vec = xx_vec(2,:);
[XX,YY] = ndgrid(xx_vec(1,:),xx_vec(2,:));
x_plot = [XX(:),YY(:)]'; % 2 x numPts^2

epsilon = zeros(size(XX));
for idx_pt = 1:numel(XX)
    x_i = x_plot(:,idx_pt);
    
    % Evaluate the expected measurement at x_i
    psi_i = triang.measurement(x_sensor,x_i);
    
    % Compute the measurement error
    err = psi-psi_i;
    
    % Evaluate the scaled log likelihood
    if do_decomp
        epsilon(idx_pt) = err'/C_d*err;
    else
        epsilon(idx_pt) = err'*C_inv*err;
    end
end
