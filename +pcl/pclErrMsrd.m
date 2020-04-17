function [epsilon,x_vec,y_vec] = pclErrMsrd(xtx,xrx,C,br_msrd,x_max,numPts,alt)
% Construct a 2-D field from -x_max to +x_max, using numPts in each
% dimension.  For each point, compute the difference in bistatic range
% between that point and the true position, normalized by the sensor
% range measurement covariance matrix C.
%
% INPUTS:
%   xtx         nDim x N transmitter position
%   xrx         nDim x M receiver positions
%   C           NM x NM covariance error matrix for bistatic range
%               measurements
%   br_msrd     NM x nMsrd matrix of measured bistatic ranges.
%   x0          nDim x 1 matrix of true emitter position
%   x_max       nDim x 1 (or scalar) vector of maximum offset from origin
%               for plotting
%   numPts      Number of test points along each dimension
%
% OUTPUTS:
%   epsilon     2-D plot of bistatic range error
%   x_vec       Vector of x coordinates for 2-D plot
%   y_vec       Vector of y coordinates for 2-D plot
%
% Nicholas O'Donoughue
% 18 Feb 2020

%% Compute the True FDOA measurement
N = size(xtx,2);
M = size(xrx,2);
R_wrapper = @(x) reshape(pcl.measurement(xtx,xrx,x),N*M,[]);

% Process measured bistatic range
% br_msrd can be specified as an NxM array, for N sensors and M range
% measurements per sensor, or an Nx1 cell array.  If the former,
% use NaN to fill in any rows with fewer range difference measurements.
if iscell(br_msrd)
    n_msmt = cellfun(@(x) numel(x),br_msrd);
    tmp = NaN*ones(numel(br_msrd),max(n_msmt));
    for i=1:numel(br_msrd)
        tmp(i,1:n_msmt(i)) = br_msrd{i};
    end
    br_msrd = tmp;
end


% Preprocess covariance matrix inverses
C_inv = inv(C);

%% Set up test points
xx_vec = x_max(:).*repmat(linspace(-1,1,numPts),2,1);
x_vec = xx_vec(1,:);
y_vec = xx_vec(2,:);
[XX,YY] = ndgrid(xx_vec(1,:),xx_vec(2,:));
x_plot = [XX(:),YY(:)]'; % 2 x numPts^2

if size(x0,1)==3
    % Add altitude
    if nargin < 6 || isempty(alt)
        alt = median(x0(3,:));
    end
    x_plot = cat(1,x_plot,alt*ones(1,size(x_plot,2)));
end

epsilon = zeros(size(XX));
for idx_pt = 1:numel(XX)
    x_i = x_plot(:,idx_pt);
    r_i = R_wrapper(x_i);
    full_err = abs(r_i-br_msrd); % N x Nmsrd, bistatic range error for each measured bistatic range at each receiver
    err = min(full_err,[],2); % For each receiver, use the closes bistatic range ellipse (smallest error)
    
%     epsilon(idx_pt) = (err'/C_d*err);
    % Compute scaled log-likelihood score
    epsilon(idx_pt) = err'*C_inv*err;
end
