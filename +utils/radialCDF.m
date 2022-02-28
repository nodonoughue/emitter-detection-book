function z = radialCDF(cov, rad, n_int_pts)
% z = radialCDF(cov, rad, n_int_pts)
%
% Compute the CDF of the radius of an elliptical bivariate normal
% distribution defined by the 2x2 matrix cov at the radius rad. The
% distribution is assumed to have zero mean.
%
% In other words, this computes the probability that the radius of a random
% 2D variable is bounded by rad, given that the 2D variable has a bivariate
% normal distribution.
%
% It does not matter if the two dimensions are correlated, or not.
% Eigenvalue decomposition will be used to orthogonalize them.
%
% The CDF is computed via numerical integration, based on the derivations 
% in:
%  Victor Chew and Ray Boyce, "Distribution of radial error in the
%  bivariate elliptical normal distribution," Technometrics, Feb 1962.
%  url: https://www.jstor.org/stable/pdf/1266181.pdf
%
% INPUTS:
%   cov     2x2 or 2x2xM covariance matrix (if 3D, the process is repeated
%           across the third dimension).
%   rad     Radius at which to compute CDF.  If multiple radii are
%           provided, then the process is repeated for each.
%   n_int_pts  (Optional) number of points to use in numerical integration
%
% OUTPUTS:
%   z       Output CDF. Dimensions are M x N, where M is the number of
%           covariance matrices provided and N is the number of radius
%           inputs.
%
% Nicholas O'Donoughue
% 3 January 2022

%% Parse Inputs

% Number of numerical integration points
if nargin < 3 || isempty(n_int_pts)
    n_int_pts = 1000;
end

[n_dim1, n_dim2, n_matrices] = size(cov);
n_radii = numel(rad);

% assert(n_dim1 == 2 && n_dim2 == 2, 'First two dimensions of cov must have size 2.');
assert(n_dim1 == n_dim2, 'First two dimensions of cov must match.');

%% Process Covariance Matrices
var_x = zeros(n_matrices, 1);
var_y = zeros(n_matrices, 1);

for idx_m = 1:n_matrices
    this_cov = cov(:,:,idx_m);

    if this_cov(1,2) == 0 && this_cov(2,1) == 0
        % Already diagonalized
        var_x(idx_m) = min(this_cov(1,1),this_cov(2,2));
        var_y(idx_m) = max(this_cov(1,1),this_cov(2,2));
    else
        % Use eigendecomposition
        lam = eig(this_cov);

        % Sort the eigenvalues
        lamSort = sort(lam,'descend');

        var_x(idx_m) = lamSort(2);
        var_y(idx_m) = lamSort(1);
    end
end

sigma_x = sqrt(var_x);
sigma_y = sqrt(var_y);

%% Check for Unstable Conditions
ratio = var_y./var_x;

mask_1drows = ratio > 1e10;
mask_infrows = isinf(ratio) || isnan(ratio);
mask_validrows = ~mask_1drows & ~mask_infrows;

% num_validrows = sum(mask_validrows);

%% Precompute Constants and Define PDF
a = (var_y + var_x)./(4 * var_x .* var_y);
b = (var_y - var_x)./(4 * var_x .* var_y);

% Restrict the a/b columns to rows that required the full calculation
a = a(mask_validrows);
b = b(mask_validrows);

% This is a slightly different form that was defined in Chew and Boyce,
% because numerical computation of the Modified Bessel Function runs into 
% precision issues as r grows large (it reaches INF very quickly).  Using
% the scaled modified bessel function of the first kind avoids this issue.
fr = @(r) r .* exp((b-a).*r.^2) .* besseli(0, b.*r.^2, 1) ./ (sigma_x.*sigma_y);

%% Set up Numerical Integration
% TODO: Define non-uniform integration interval based on shape of fr
%
% For example, if rad > 5 * cep50(cov), then CDF = .999 and there is very
% little contribution.
% Also, if sigma_x / sigma_y >> 1, we can start to approximate as a 1D
% problem.
dr = max(rad)/(n_int_pts-1);
r_vec = 0:dr:max(rad);

% Be careful to ensure that r_vec is a row vector, since a and b are column
% vectors (one entry for each of the provided covariance matrices).
pdf = fr(r_vec); % num_validrows x n_int_points

z = zeros(n_matrices, n_radii);
for idx_r = 1:n_radii
    out_idx = find(r_vec > rad(idx_r),1,'first') - 1;
    z(valid_mask, idx_r) = sum(pdf(:,1:out_idx),2)*dr;
end

%% Handle the invalid rows
z(mask_1drows,:) =  utils.normcdf(rad,0,sigma_y(mask_1drows)) ...
                   -utils.normcdf(-rad,0,sigma_y(mask_1drows));
z(mask_infrows,:) = 0;
