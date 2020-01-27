function t = backtrackingLineSearch(f,x,grad,del_x,alpha,beta)
% t = backtrackingLineSearch(f,x,grad,del_x,alpha,beta)
%
% Performs backtracking line search according to algorithm 9.2 of
% Stephen Boyd's, Convex Optimization
%
% Inputs:
%
%   f       Function handle to minimize
%   x       Current estimate of x
%   grad    Gradient of f at x
%   del_x   Descent direction
%   alpha   Constant between 0 and 0.5
%   beta    Constant between 0 and 1
%
% Outputs:
%
%   t       Optimal step size for the current estimate x.
%
% Nicholas O'Donoughue
% 1 July 2019

% Initialize the search parameters and direction
t = 1;
startingVal = f(x);
slope = grad'*del_x;

% Make sure the starting value is large enough
while f(x+t*del_x) < startingVal+alpha*t*slope
    t = 2*t;
end

% Conduct the backtracking line search
while f(x+t*del_x) > startingVal+alpha*t*slope
    t = beta*t;
end
