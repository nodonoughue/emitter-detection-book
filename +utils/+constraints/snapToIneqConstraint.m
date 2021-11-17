function x_valid = snapToIneqConstraint(x, b)
% function x_valid = snapToIneqConstraint(x, b)
%
% Apply the inequality constraints in the function handle b to the position 
% x.
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
% and 5.11, and Figure 5.6.  This is intended for use in iterative solvers.
%
% INPUTS:
%   x       nDim x N set of input variables against which to apply the
%           equality constraints.
%   b       nConst x 1 array of function handles representing equality
%           constraints.
%
% OUTPUTS:
%   x_valid nDim x N matrix of validated input variables
%
% Nicholas O'Donoughue
% 14 November 2021


%% Parse Inputs
[~, N] = size(x); % Specifying 2 outputs collapses any higher 
                     % dimensions of x
nConst = numel(b);

%% Multiple Variables Provided; Use Recursion to solve
if N > 1
    % Use arrayfun across the 2nd dimension of x. The response is an N x 1
    % cell array, each cell contains an nDim x 1 vector.
    x_cell = arrayfun(@(n) snapToIneqConstraint(x(:,n), b), ...
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
    this_b = b{idx};
    
    % Test the current value (what is the error, and what scale parameter
    % are needed to force equality)
    [epsilon, this_x_valid] = this_b(x_valid);
    
    if epsilon > 0
        % Inequality constraint is broken, apply scale factor
        x_valid = this_x_valid;
    end
end
