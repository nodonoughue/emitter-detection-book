function ell = drawEllipseFromFoci(xf,xpt,numPts,dist_offset)
% Draw an ellipse from the foci in xf (nDim x 2) and a point on the
% ellipse xpt (nDim x 1).  
%
% If xpt has more than one column, a best fit ellipse will be drawn based
% on the average sum distance between the points provided and the foci.
%
% There is an option third argument, numPts, that specifies how many point
% to draw in the ellipse (evenly spread in angle about the center point).
%
% There is an optional fourth argument, dist_offset, that will provide a
% delta to the sum distance computed, to allow for the ellipse to be drawn
% smaller or larger.
%
% Nicholas O'Donoughue
% 24 February 2020

% Default density of points
if nargin < 3 || isempty(numPts)
    numPts = 1000;
end

if nargin < 4 || isempty(dist_offset)
    dist_offset = 0;
end

% Find the rotation angle from the x-axis
rotAngle = atan2(diff(xf(2,:)),diff(xf(1,:)));
foci_dist = sqrt(sum(abs(xf(:,1)-xf(:,2)).^2));
x_ctr = mean(xf,2); % Center point

% Find the major axis length, based on distance from the point in xpt to
% the two foci
sum_dist = sqrt(sum(abs(xpt-xf(:,1)).^2,1))+sqrt(sum(abs(xpt-xf(:,2)).^2,1))+dist_offset;

% Check if multiple points provided
if numel(sum_dist)>1
    % Take the average distance
    sum_dist = mean(sum_dist);
end

% Find the major and minor axis lengths
major_axis = sum_dist/2; % distance from the center point to each vertex on major axis
minor_axis = sqrt(major_axis^2 - (foci_dist/2)^2);

% Define Error Ellipse in rotated coordinate frame
th = linspace(0,2*pi,numPts); % Angle from center point, in rotated frame
ellipseX = major_axis*cos(th);         % Major axis coordinates
ellipseY = minor_axis*sin(th);         % Minor axis coordinates
    
% Rotate the ellipse to standard reference frame
R = [cos(rotAngle), -sin(rotAngle); sin(rotAngle), cos(rotAngle)];
xx = R*[ellipseX;ellipseY]; % Store as a 2 x N matrix of positions

% Apply bias; so the ellipse is centered on the input x
ell = x_ctr+xx;


