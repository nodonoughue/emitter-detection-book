function err_proj = projectError(err, x_pt, x_obs)
% err_proj = projectError(err, x_pt, x_obs)
%
% Project a 3D error (err) covariance matrix perpendicular to the line of
% sight from x_obs to x_pt.
%
% Error covariance is either 3 x 3 (x, y, z) or 6 x 6 (x, y, z, vx, vy, vz)
%
% INPUTS:
%   err         3 x 3 x N datacube of covariance matrices for 
%               N different test cases.
%   x_pt        3 x 1 or 3 x N position vector for the target
%   x_obs       3 x 1 or 3 x N position vector for the observer
%
% OUTPUTS:
%   err_proj    3 x 3 x N datacube of projected covariance
%               matrices for N different test cases.
%
% Nicholas O'Donoughue
% 8 Nov 2021

%% Parse Inputs
[n_dim, n_case_pt] = size(x_pt);
[n_dim_obs,n_case_obs] = size(x_obs);
[n_dim_cov,n_dim2_cov,n_case] = size(err);

assert(n_dim_cov == n_dim2_cov, 'Error covariance mismatch.');
assert((n_case_pt == 1 || n_case_pt == n_case) && ...
       (n_case_obs == 1 || n_case_obs == n_case), 'Unequal number of test cases.');
   
assert(n_dim == n_dim_obs,'Error: target and observation vector must have matching size.');
assert(n_dim_cov == n_dim, 'Error: covariance matrix dimensions must match target position.');


%% Define LOS Vector and Projection Matrix
dx = x_pt - x_obs;  % 3 x 1 or 3 x N
u  = dx / norm(dx); % unit vector

proj = u*u';
proj_ortho = eye(n_dim) - proj;

%% Project Errors
err_cell = arrayfun(@(n) proj_ortho*squeeze(err(1:n_dim,1:n_dim,n))*proj_ortho, 1:n_case,'UniformOutput',false);
    % Project each of the N matrices orthogonally to u
    % Result is a cell-array of n_dim x n_dim matrices.

err_proj = cell2mat(reshape(err_cell,1,1,n_case));