function x_valid = snapToEqConstraint(x, a, tol)
% function x_valid = snapToEqConstraint(x, a, tol)
%
% Apply the equality constraints in the function handle a to the position 
% x, subject to a tolerance tol.
%
% If multiple constraints are provided, then they are applied in order. In
% this manner, only the final one is guaranteed to be satisfied by the
% output x_valid, as application of later constraints may violate earlier
% ones.
%
% If the variable x has non-singleton dimensions (beyond the first), then
% this operation is repeated across them; the first dimensions is assumed
% to be used for vector inputs x (such as a 2D or 3D target position).
%
% The projection operation is carried out as discussed in equations 5.10
% and 5.11.  This is intended for use in iterative solvers.
%
% INPUTS:
%   x       nDim x N set of input variables against which to apply the
%           equality constraints.
%   a       nConst x 1 array of function handles representing equality
%           constraints.
%   tol     Tolerance of the equality constraints.  If left unspecified,
%           or empty, then the default value is 1e-6.
%
% OUTPUTS:
%   x_valid nDim x N matrix of validated input variables
%
% Nicholas O'Donoughue
% 14 November 2021

%% Parse Inputs
[~, N] = size(x); % Specifying 2 outputs collapses any higher 
                     % dimensions of x
if isa(a,'function_handle')
    % Wrap it in a cell array for indexing
    a = {a};
end
nConst = numel(a);

if nargin < 3 || ~exist('tol','var') || isempty(tol)
    tol = 1e-6;
end

%% Multiple Variables Provided; Use Recursion to solve
if N > 1
    % Use arrayfun across the 2nd dimension of x. The response is an N x 1
    % cell array, each cell contains an nDim x 1 vector.
    x_cell = arrayfun(@(n) snapToEqConstraint(x(:,n), a, tol), ...
                      'UniformOutput',false); 

    % Reshape the cell array, then use cell2mat.  Result is an nDim x N
    % matrix.  Use another reshape command to ensure that the shape of any
    % higher dimensions matches the input x.
    x_valid = reshape(cell2mat(reshape(x_cell,1,N)), size(x));
    return;
end

%% Single Position Variable
x_valid = x; % Initialize the output
for idx=1:nConst
    % Isolate the current equality constraint (function handle)
    this_a = a{idx};
    if isempty(this_a), continue, end
    if ~isa(this_a,'function_handle'), continue, end

    % Test the current value (what is the error, and what scale parameter
    % are needed to force equality)
    [epsilon, this_x_valid] = this_a(x_valid);

    if abs(epsilon) > tol
        % Equality constraint is broken, apply scale factor
        x_valid = this_x_valid;
    end

    % Note that we're updating x_valid to match each constraint
    % one-by-one.  It is likely, if there are multiple constraints,
    % that when we align to the last, we are no longer aligned to the
    % first.
end
