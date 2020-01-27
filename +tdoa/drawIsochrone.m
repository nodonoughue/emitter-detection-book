function iso = drawIsochrone(x1,x2,rdiff,numPts,maxOrtho)
% iso = drawIsochrone(x1,x2,rdiff,numPts,maxOrtho)
%
% Finds the isochrone with the stated range difference from points x1
% and x2.  Generates an arc with 2*numPts-1 points, that spans up to
% maxOrtho distance from the intersection line of x1 and x2
%
% Inputs:
%   x1          Sensor position 1
%   x2          Sensor position 2
%   rdiff       Desired isochrone range difference
%   numPts      Number of points to draw
%   maxOrtho    Maximum offset (orthogonal to the line connecting x1 and
%               x2) at which points should be drawn
%
% Outputs:
%   iso         2 x numPts array of isochrone coordinates
%
% Nicholas O'Donoughue
% 1 July 2019

% Generate pointing vectors u and v in rotated coordinate space
%  u = unit vector from x1 to x2
%  v = unit vector orthogonal to u
rotMat = [0,1; -1,0];
R = utils.rng(x1,x2);
u = (x2-x1)/R;
v = rotMat*u;

xProj = [u v];

% Position of reference points in uv-space
x1uv = [0;0];
x2uv = [R;0];

% Initialize isochrone positions in uv-space
vv = linspace(0,maxOrtho,numPts);
uu = zeros(size(vv));
uu(1) = (R-rdiff)/2;
xuv = [uu; vv];

% Integrate over points, search for isochrone position
maxIter = 10000;

for i=1:numPts
    if i>1
        xuv(1,i) = xuv(1,i-1); % Initialize u position with previous value
    end
    
    offset = R;
    numIter = 1;
    
    while offset > R/10000 && numIter <= maxIter
        numIter = numIter + 1;
        
        rdiff0 = utils.rngDiff(xuv(:,i),x1uv,x2uv);
        offset = rdiff0-rdiff;
        
        xuv(1,i) = xuv(1,i)+offset;
    end
end

% Isochrone is symmetric about u axis
    % Flip the u axis, flip and negate v
xuv = [fliplr([xuv(1,2:end);-xuv(2,2:end)]),xuv];

% Convert to x/y space and re-center at origin
iso = xProj*xuv + x1;
    