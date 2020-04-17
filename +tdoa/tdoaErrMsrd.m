function [epsilon,x_vec,y_vec] = tdoaErrMsrd(xs,C,r_msrd,x_max,numPts,alt)
% [epsilon,x_vec,y_vec] = tdoaErrMsrd(xs,C,r_msrd,x_max,numPts,alt)
%
% Construct a 2-D field from -x_max to +x_max, using numPts in each
% dimension.  For each point, compute the TDOA solution for each sensor
% against the reference (the first sensor), and compare to the TDOA
% solution from the true emitter position.
%
% INPUTS:
%   xs          nDim x N matrix of sensor positions
%   C           N x N covariance error matrix for Range Difference
%               measurements
%   x0          nDim x 1 matrix of true emitter position
%   x_max       nDim x 1 (or scalar) vector of maximum offset from origin
%               for plotting
%   numPts      Number of test points along each dimension
%
% OUTPUTS:
%   epsilon     2-D plot of TDOA error
%   x_vec       Vector of x coordinates for 2-D plot
%   y_vec       Vector of y coordinates for 2-D plot
%
% Nicholas O'Donoughue
% 18 Feb 2020

%% Compute the True FDOA measurement
xx0 = xs(:,1:end-1); % test sensors
xxN = xs(:,end); % reference sensor
R0 = @(x) sqrt(sum(abs(x-xx0).^2,1)); % eq 11.4
RN = @(x) sqrt(sum(abs(x-xxN).^2,1)); % eq 11.5
dr  = @(x) reshape(R0(x)-RN(x),size(xx0,2),[]);

% Process measured range difference
% r_msrd can be specified as an NxM array, for N sensor pairs and M range
% difference measurements per pair, or an Nx1 cell array.  If the former,
% use NaN to fill in any rows with fewer range difference measurements.
if iscell(r_msrd)
    n_msmt = cellfun(@(x) numel(x),r_msrd);
    tmp = NaN*ones(numel(r_msrd),max(n_msmt));
    for i=1:numel(r_msrd)
        tmp(i,1:n_msmt(i)) = r_msrd{i};
    end
    r_msrd = tmp;
end

% Preprocess covariance matrix inverses
C_out = C(1:end-1,1:end-1) + C(end,end);
% C_d = decomposition(C_out);
C_inv = inv(C_out);

%% Set up test points
xx_vec = x_max(:).*repmat(linspace(-1,1,numPts),2,1);
x_vec = xx_vec(1,:);
y_vec = xx_vec(2,:);
[XX,YY] = ndgrid(xx_vec(1,:),xx_vec(2,:));
x_plot = [XX(:),YY(:)]'; % 2 x numPts^2

if size(xs,1)==3
    % Add altitude
    if nargin < 6 || isempty(alt)
        alt = median(xs(3,:));
    end
    x_plot = cat(1,x_plot,alt*ones(1,size(x_plot,2)));
end

epsilon = zeros(size(XX));
for idx_pt = 1:numel(XX)
    x_i = x_plot(:,idx_pt);
    r_i = dr(x_i);
    err = min(abs(r_msrd-r_i),[],2); % Closest RDOA line
%     epsilon(idx_pt) = (err'/C_d*err);
    epsilon(idx_pt) = err'*C_inv*err;
end
