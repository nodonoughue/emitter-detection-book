function ell = drawErrorEllipse(x,C,numPts,confInterval)
% ell = drawErrorEllipse(x,C,numPts,confInterval)
%
% Compute and return the error ellipse coordinates with numPts 
% samples around the ellipse, for an error centered at the location x and
% with covariance matrix C.  The confidence interval specifies what
% percentage of random errors will fall within the error ellipse.
%
% Inputs:
%   x           Center position of ellipse (2-dimensional coordinates)
%   C           Error covariance matrix (2x2)
%   numPts      Number of points to use to draw the ellipse [default = 100]
%   confInterval    Desired confidence interval, must be one of:
%                       1  -- 1 standard deviation (67%) [default]
%                       50 -- 50 percent
%                       90 -- 90 percent
%                       95 -- 95 percent
%
% Outputs:
%   ell         2 x N vector of coordinates defining the error ellipse
%
% Nicholas O'Donoughue
% 1 July 2019

if nargin < 4 || ~exist('confInterval','var') || isempty(confInterval)
    confInterval = 1; % 1 std. dev is the default
end

if nargin < 3 || ~exist('numPts','var') || isempty(numPts)
    numPts = 100;
end

% Eigenvector analysis to identify major/minor axes rotation and length
[V,Lam] = eig(C);
lam = diag(Lam); % Pull eigenvalue vector from diagonal matrix Lam

% Sort the eigenvalues
[lamSort,iSort] = sort(lam,'descend'); 

% Major Axis
lamMax = lamSort(1);
vMax = V(:,iSort(1)); 

% Minor Axis
lamMin = lamSort(end);
% vMin = V(:,iSort(end)); % Not actually used -- but left in for
                          % completeness

% Compute the rotation angle
rotAngle = atan2(vMax(2),vMax(1));
rotAngle = mod(rotAngle+pi,2*pi)-pi; % ensure angle is on interval [-pi,pi]

% Lookup scale factor from confidence interval
if confInterval==1
    gamma = 1; % 1 sigma
elseif confInterval==50
    gamma = 1.386;
elseif confInterval==90
    gamma=4.601;
elseif confInterval==95
    gamma = 5.991; % 95% CI
else
    gamma = 1;
    fprintf('Confidence Interval not recognized, using 1 standard deviation...\n');
end

% Define Error Ellipse in rotated coordinate frame
th = linspace(0,2*pi,numPts); % Angle from center point, in rotated frame
a = sqrt(gamma*lamMax);       % Major axis length
b = sqrt(gamma*lamMin);       % Minor axis length
    
ellipseX = a*cos(th);         % Major axis coordinates
ellipseY = b*sin(th);         % Minor axis coordinates
    
% Rotate the ellipse to standard reference frame
R = [cos(rotAngle), -sin(rotAngle); sin(rotAngle), cos(rotAngle)];
xx = R*[ellipseX;ellipseY]; % Store as a 2 x N matrix of positions

% Apply bias; so the ellipse is centered on the input x
ell = x+xx;