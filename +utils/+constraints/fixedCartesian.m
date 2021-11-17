function [a, a_grad] = fixedCartesian(type, varargin)
% [a, a_grad] = fixedCartesian(type, bound)
% [a, a_grad] = fixedCartesian(type, x0, u)
%
% Generate a set of constraints and constraint gradient functions for fixed
% solution constraints on a cartesian grid.
%
% If type is 'x', 'y', or 'z', creates a fixed bound along one of those
% axes, according to the scalar input bound.
%
% If type is 'linear', creates a linear bound defined by a point x0 and
% pointing vector u.
%
% For example:
%   utils.constraints.fixedCartesian('x',5)
% will generate a constraint that fixes the first dimension to 5, but
% allows the second and (if defined) third dimension to vary unconstrained.
%
%
% Inputs:
%   type    String ('x', 'y', 'z', or 'linear')
%   bound   Fixed bound for x, y, or z type constraints
%   x0      Fixed point for linear type constraints
%   u       Pointing vector for linear type constraints
%
% Outputs:
%   a       Array of constraint function handles, each returning an epsilon
%           (or distance from equality), and a scale term (against which 
%           the input x should be multiplied to achieve equality).
%   a_grad  Array of constraint gradient function handles
%
% Nicholas O'Donoughue
% 14 November 2021

%% Parse Inputs
assert(nargin >= 2, 'Not enough inputs');

switch lower(type)
    case 'x'
        a = @(x) fixedCartX(x, varargin{1});
        a_grad = @(x) [1;zeros(size(x,1)-1,1)]*ones(1,size(x,2));
    case 'y'
        a = @(x) fixedCartY(x, varargin{1});
        a_grad = @(x) [0;1;zeros(size(x,1)-2,1)]'*ones(1,size(x,2));
    case 'z'
        a = @(x) fixedCartZ(x, varargin{1});
        a_grad = @(x) [0;0;1]*ones(1,size(x,2));
    case 'linear'
        assert(nargin > 2, 'Not enough inputs.');
        a = @(x) fixedCartLinear(x, varargin{1}, varargin{2});
        a_grad = fixedCartLinearGrad(x, varargin{1}, varargin{2});
    otherwise
        error('Unrecognized constraint type. Please use ''x'', ''y'', ''z'', or ''linear''');
end



%% Subordinate Functions
    function [epsilon, x_valid] = fixedCartX(x, x_c)
        epsilon = x(1,:) - x_c;

        x_valid = x;
        x_valid(1,:) = x_c;
    end

    function [epsilon, x_valid] = fixedCartY(x, y_c)
        epsilon = x(2,:) - y_c;

        x_valid = x;
        x_valid(2,:) = y_c;
    end

    function [epsilon, x_valid] = fixedCartZ(x, z_c)
        epsilon = x(3,:) - z_c;

        x_valid = x;
        x_valid(3,:) = z_c;
    end

    function [epsilon, x_valid] = fixedCartLinear(x, x_0, u)
        u = u/norm(u); % Make sure u is unit norm
        n_dim = size(u,1);
        P = eye(n_dim) - u(:)*u(:)'; % projection, orthogonal to u

        % Component of x-x_0 that is perpendicular to u
        y = P*(x-x_0);

        % Distance of x from the vector x0+delta*u
        epsilon = sqrt(y'*y);

        % To project x onto the constraint, simply subtract y, the
        % component of the vector x-x0 that is perpendicular to it
        x_valid = x - y;
    end

    function grad = fixedCartLinearGrad(x, x_0, u)
        u = u/norm(u); % Make sure u is unit norm
        n_dim = size(u,1);
        P = eye(n_dim) - u(:)*u(:)'; % projection, orthogonal to u

        grad = 2*P*(x-x_0);
    end
end