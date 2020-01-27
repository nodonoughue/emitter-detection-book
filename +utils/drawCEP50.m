function x0 = drawCEP50(x,C,numPts)
% x0 = drawCEP50(x,C,numPts)
%
% Return the (x,y) coordinates of a circle with radius given by the CEP50
% of the covariance matrix C, centered on the point x, with numPts points
% around the circle
%
% Inputs:
%   x           Center position for CEP circle (2-dimensional)
%   C           Covariance matrix (2x2)
%   numPts      Number of points to use to draw circle
%
% Outputs:
%   x0          2xnumPts array of coordiantes defining CEP circle
%
% Nicholas O'Donoughue
% 1 July 2019

% Find the CEP; the radius of the circle to draw
cep = utils.computeCEP50(C);

% Initialize the angular indices of the circle
th = linspace(0,2*pi,numPts);

% Convert from r,th to cartesian coordinates
xx = cep*cos(th);
yy = cep*sin(th);

% Offset the cartesian coordinates
x0 = bsxfun(@plus,x,[xx;yy]);  
