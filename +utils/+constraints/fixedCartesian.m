function [a, a_grad] = fixedCartesian(x_c, y_c, z_c)
% [a, a_grad] = fixedCartesian(x_c, y_c, z_c)
%
% Generate a set of constraints and constraint gradient functions for fixed
% solution constraints on a cartesian grid.  Inputs are the fixed
% constraints in x, y, and z, respectively.
%
% To allow a dimension to vary unconstrained, leave the corresponding input
% blank.
%
% For example:
%   utils.constraints.fixedCartesian(5,[],[])
% will generate a constraint that fixes the first dimension to 5, but
% allows the second and (if defined) third dimension to vary unconstrained.
%
% Similarly,
%   utils.constraints.fixedCartesian([],3,[])
% will fix the second dimension, while allowing the first and (if defined)
% third to vary.
%
% To constrain two dimensions, supply two variables, such as:
%   utils.constraints.fixedCartesian(1,[],6)
% or
%   utils.constraints.fixedCartesian([],2,4)
%
% Inputs:
%   x_c     Constraint for the first dimension (empty if unconstrained)
%   y_c     Constraint for the second dimension (empty if unconstrained)
%   z_c     Constraint for the third dimension (empty if unconstrained)
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
do_x = nargin > 0 && exist('x_c','var') && ~isempty(x_c);
do_y = nargin > 1 && exist('y_c','var') && ~isempty(y_c);
do_z = nargin > 2 && exist('z_c','var') && ~isempty(z_c);

assert(~(do_x & do_y & do_z), 'At least one dimension must be left to vary.')

%% Fixed Constraint on z
if do_x
    a_x = @(x) fixedCartX(x, x_c);
    a_grad_x = @(x) [1,0,0]'*ones(1,size(x,2));
else
    a_x = [];
    a_grad_x = [];
end

%% Fixed Constraint on y
if do_y
    a_y = @(x) fixedCartY(x, y_c);
    a_grad_y = @(x) [0,1,0]'*ones(1,size(x,2));
else
    a_y = [];
    a_grad_y = [];
end

%% Fixed Constraint on z
if do_z
    a_z = @(x) fixedCartZ(x, z_c);
    a_grad_z = @(x) [0,0,1]'*ones(1,size(x,2));
else
    a_z = [];
    a_grad_z = [];
end

%% Compile outputs
a = {a_x, a_y, a_z};
a_grad = {a_grad_x, a_grad_y, a_grad_z};

%% Subordinate Functions
    function [epsilon, scale] = fixedCartX(x, x_c)
        [n_dim,n_x] = size(x);
        epsilon = x(1,:) - x_c;
        scale = [x_c./x(1,:); ones(n_dim-1,n_x)];
    end

    function [epsilon, scale] = fixedCartY(x, y_c)
        [n_dim,n_x] = size(x);
        epsilon = x(2,:) - y_c;
        scale = [ones(1,n_x); y_c./x(2,:); ones(n_dim-2,n_x)];
    end

    function [epsilon, scale] = fixedCartZ(x, z_c)
        [~,n_x] = size(x);
        epsilon = x(3,:) - z_c;
        scale = [ones(2,n_x); z_c./x(3,:)];
    end
end