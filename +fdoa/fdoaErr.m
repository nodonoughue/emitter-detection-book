function [epsilon,x_vec,y_vec] = fdoaErr(xs,xs_dot,C,x0,x_max,numPts)
% [epsilon,x_vec,y_vec] = fdoaErr(xs,xs_dot,C,x0,x_max,numPts)
%
% Construct a 2-D field from -x_max to +x_max, using numPts in each
% dimension.  For each point, compute the FDOA solution for each sensor
% against the reference (the first sensor), and compare to the FDOA
% solution from the true emitter position.
%
% INPUTS:
%   xs          nDim x N matrix of sensor positions
%   xs_dot      nDim x N matrix of sensor velocities
%   x0          nDim x 1 matrix of true emitter position
%   x_max       nDim x 1 (or scalar) vector of maximum offset from origin
%               for plotting
%   numPts      Number of test points along each dimension
%
% OUTPUTS:
%   epsilon     2-D plot of FDOA error
%   x_vec       Vector of x coordinates for 2-D plot
%   y_vec       Vector of y coordinates for 2-D plot
%
% Nicholas O'Donoughue
% 1 July 2019

%% Compute the True FDOA measurement
xx0 = xs(:,1:end-1); % test sensors
xxN = xs(:,end); % reference sensor
vv0 = xs_dot(:,1:end-1);
vvN = xs_dot(:,end);
R0 = @(x) sqrt(sum(abs(x-xx0).^2,1)); % eq 11.4
RN = @(x) sqrt(sum(abs(x-xxN).^2,1)); % eq 11.5
R0_dot = @(x) sum(vv0.*(x-xx0),1)./R0(x);   % eq 12.3
RN_dot = @(x) sum(vvN.*(x-xxN),1)./RN(x);   % eq 12.3
r_dot  = @(x) R0_dot(x)'-RN_dot(x)';  % eq 13.5

r = r_dot(x0); % True Range Difference

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
    r_i = r_dot(x_i);
    err = r-r_i;
    epsilon(idx_pt) = (err'/C_d*err);
end
