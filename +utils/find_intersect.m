function x = find_intersect(x0,psi0,x1,psi1)
% x = find_intersect(x0,psi0,x1,psi1)
%
% Find the intersection of two lines of bearing with origins at x0 and x1,
% and bearing given by psi0 and psi1.
%
% INPUTS:
%
%   x0          2-D position of the first sensor
%   psi0        Bearing of the first line, which begins at x0
%   x1          2-D position of the second sensor
%   psi1        Bearing of the second line, which begins at x0
%
% OUTPUTS:
%
%   x           2-D position of the intersection point
%
% Nicholas O'Donoughue
% 1 July 2019

%% Fine slope and intercept for each line of bearing
m0 = sin(psi0)/cos(psi0);
m1 = sin(psi1)/cos(psi1);

b0 = x0(2) - m0*x0(1);
b1 = x1(2) - m1*x1(1);

if (isinf(m0) && isinf(m1)) || abs(m0-m1) < 1e-12
    % They're parallel, effectively
    fprintf('Slopes are almost parallel; the solution is ill-conditioned.');
    x = x0;
    return;
end

%% Check Boundary Cases
if abs(cos(psi0))<1e-10
    % First LOB is parallel to y axis; x is fixed
    x(1) = x0(1);

    % Use slope/intercept definition of second LOB to solve for y
    x(2) = m1 * x(1) + b1;
elseif abs(cos(psi1))<1e-10
    % Same issue, but with the second LOB being parallel to y axis
    x(1) = x1(1);
    x(2) = m0 * x(1) + b0;
else
    % Find the point of intersection
    x(1) = (b0-b1)/(m1-m0);
    x(2) = m1*x(1) + b1;
end
    
