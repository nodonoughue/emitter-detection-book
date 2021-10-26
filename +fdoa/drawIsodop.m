function xy_iso = drawIsodop(x1,v1,x2,v2,vdiff,numPts,maxOrtho)
% xy_iso = drawIsodop(x1,v1,x2,v2,vdiff,numPts,maxOrtho)
%
% Finds the isochrone with the stated range rate difference from points x1
% and x2.  Generates an arc with 2*numPts-1 points, that spans up to
% maxOrtho distance from the intersection line of x1 and x2
%
% Inputs:
%   x1          Position of first sensor (Ndim x 1) [m]
%   v1          Velocity vector of first sensor (Ndim x 1) [m/s]
%   x2          Position of second sensor (Ndim x 1) [m]
%   v2          Velocity vector of second sensor (Ndim x 1) [m/s]
%   vdiff       Desired velocity difference [m/s]
%   numPts      Number of points to compute
%   maxOrtho    Maximum offset from line of sight between x1 and x2 [m]
%
% Outputs:
%   x_iso       2 x numPts iso doppler contour [m]
%
% Nicholas O'Donoughue
% 1 July 2019

% Set frequency to 3e8, so that c/f_0 is unity, and output of utils.dopDiff
% is velocity difference [m/s]
f_0 = 3e8;

%% Set up test points
xx_vec = maxOrtho(:).*repmat(linspace(-1,1,numPts),2,1);
x_vec = xx_vec(1,:);
y_vec = xx_vec(2,:);
[XX,YY] = meshgrid(xx_vec(1,:),xx_vec(2,:));
x_plot = [XX(:),YY(:)]'; % numPts^2 x 2

df_plot = utils.dopDiff(x_plot,[0 0]',x1,v1,x2,v2,f_0);
    % numPts^2 x (N-1)
    
    
%% Compute contour
fig00=figure;
[cc,hh] = contour(x_vec,y_vec,reshape(df_plot,numPts,numPts),[vdiff vdiff]);
% Close the figure generated
close(fig00);

% The contour output is formatted as a 2xN matrix with contours repeated
% in the column dimension.  Each contour starts with a header column:
%   column 1, row 1 -- contour level
%   column 1, row 2 -- number of points
%   column 2:N+1, row 1 -- x points
%   column 2:N+1, row 2 -- y points
% This is repeated for as many contours as are defined

xy_iso = [];

while numel(cc) > 0
    % level = cc(1,1);
    n_points = cc(2,1);

    xy_iso=cat(2,xy_iso,cc(:,1+(1:n_points)),[NaN;NaN]);
   
    cc = cc(:,n_points+2:end);
end

% x_iso = xy_iso(1,:);
% y_iso = xy_iso(2,:);

%% Filter points out of bounds
% out_of_bounds = abs(x_iso) > maxOrtho | abs(y_iso) > maxOrtho;
% x_iso = x_iso(~out_of_bounds);
% y_iso = y_iso(~out_of_bounds);

% xy_iso = cat(1, x_iso(:)', y_iso(:)');