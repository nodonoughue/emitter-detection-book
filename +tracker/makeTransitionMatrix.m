function F = makeTransitionMatrix(update_rate,num_states,num_dims)
% F = makeTransitionMatrix(update_rate,num_states,num_dims)
%
% State matrix is organized first with position (for the number of
% dimensions specified), then with velocity and acceleration, also with the
% number of dimensions specified.
%
% A 3-state, 3D state vector is:
%   [x, y, z, vx, vy, vz, ax, ay, az]
%
%
% INPUTS:
%   num_states      Number of states (per dimension)
%                   1 = position
%                   2 = position/velocity
%                   3 = position/velocity/accel
%   num_dims        Number of dimensions
%
% OUTPUT:
%



% 3x3 Transition Matrix (single-dimension)
if num_states == 3
    f_single_dim = @(t) [1, t, .5*t.^2; 0, 1, t;0, 0, 0];
elseif num_states == 2
    f_single_dim = @(t) [1, t; 0 1];
else
    f_single_dim = @(t) 1;
end

% Dimension Distribution Matrix
k = eye(num_dims);

% Use the Kronecker product to distribute the single-dimension transition
% matrices.
F = kron(k,f_single_dim(update_rate));

