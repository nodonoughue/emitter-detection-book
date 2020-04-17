function [epsilon,x_vec,y_vec,epsilon_r,epsilon_rr] = pclErr(xtx,xrx,C,x0,x_max,numPts,alt,vtx,vrx,Crr,v0)
% Construct a 2-D field from -x_max to +x_max, using numPts in each
% dimension.  For each point, compute the difference in bistatic range
% between that point and the true position, normalized by the sensor
% range measurement covariance matrix C.
%
% INPUTS:
%   xtx         nDim x N transmitter position
%   xrx         nDim x M receiver positions
%   C           NM x NM covariance error matrix for bistatic range
%               measurements (for each of the NM sensor pairs)
%               Note: measurements will vary across all transmitters for
%               the first receiver, then repeat for each receiver.  Cov
%               ariance matrix must be similarly constructed.
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

do_doppler = nargin >= 8;

%% Compute the True FDOA measurement
N = size(xtx,2);
M = size(xrx,2);
nTgt = size(x0,2);
r_true = reshape(pcl.measurement(xtx,xrx,x0),N*M,[]);
R_wrapper = @(x) reshape(pcl.measurement(xtx,xrx,x),N*M,[]) - r_true;

% Preprocess covariance matrix inverses
Cr_inv = inv(C);

if do_doppler
    % Handle defaults
    if isempty(vtx)
        vtx = zeros(size(xtx));
    end
    
    if isempty(vrx)
        vrx = zeros(size(xrx));
    end
    
    if isempty(v0)
        v0 = zeros(size(x0));
    end
    
    if isempty(Crr)
        Crr = eye(size(C,1));
    end
    
    % Range rate true
    [~,rr_true] = pcl.measurement(xtx,xrx,x0,vtx,vrx,v0);
    rr_true = reshape(rr_true,N*M,[]);
    
    % Range rate residual error wrapper
    RR_wrapper = @(x,rr_est) pcl.computeRngRateResidual(xtx,xrx,x,vtx,vrx,rr_est);
    
    % Add range rate variance
    Crr_inv = inv(Crr);
    C_inv = blkdiag(Cr_inv,Crr_inv);
    %C_inv = inv(Crr);
else
    C_inv = Cr_inv;
end

%% Set up test points
xx_vec = x_max(:).*repmat(linspace(-1,1,numPts),2,1);
x_vec = xx_vec(1,:);
y_vec = xx_vec(2,:);
[XX,YY] = ndgrid(xx_vec(1,:),xx_vec(2,:));
x_plot = [XX(:),YY(:)]'; % 2 x numPts^2

if size(x0,1)==3
    % Add altitude
    if nargin < 7 || isempty(alt)
        alt = median(x0(3,:));
    end
    x_plot = cat(1,x_plot,alt*ones(1,size(x_plot,2)));
end

epsilon = zeros(size(XX));
epsilon_r = zeros(size(XX));
epsilon_rr = zeros(size(XX));

for idx_pt = 1:numel(XX)
    x_i = x_plot(:,idx_pt);
    rng_err = R_wrapper(x_i); % NM x nTgt
    
    % Find the closest solution for each measurement
    [min_rng_err,idx_tgt] = min(abs(rng_err),[],2);
    
    if do_doppler
        % Grab the range-rate measurement from the closest range-rate
        % solutions
        rr_est = rr_true(sub2ind([N*M,nTgt],(1:(N*M))',idx_tgt));
        rr_err = RR_wrapper(x_i,rr_est); % NM x nTgt
        
        % Grab the same residuals
        err = [min_rng_err;rr_err];
    else
        err = min_rng_err;
    end
    
%     epsilon(idx_pt) = (err'/C_d*err);
    % Compute scaled log-likelihood score
    epsilon(idx_pt) = err'*C_inv*err;
    epsilon_r(idx_pt) = min_rng_err'*Cr_inv*min_rng_err;
    if do_doppler
        epsilon_rr(idx_pt) = rr_err'*Crr_inv*rr_err;
    else
        epsilon_rr(idx_pt) = 0;
    end
end
