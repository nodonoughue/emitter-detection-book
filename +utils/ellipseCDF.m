function cdf = ellipseCDF(cov,rad)
% cdf = ellipseCDF(cov,rad)
%
% Compute the CDF of an elliptical distribution defined by cov (in x and
% y coordinates), at the radius rad.
%
% In other words, this computes the probability that the radius of a random
% 2D variable is bounded by rad, given that the 2D variable has a bivariate
% normal (Gaussian) distribution with zero mean and covariance matrix cov.
%
% INPUTS:
%   cov     2x2 covariance matrix
%   rad     desired radius
%
% OUTPUTS:
%   cdf     Probability that a random draw falls within a circle of radius
%           rad
%
% Relies on the fact that the marginal distribution of the radius of a
% bivariate Gaussian with no correlation between x and y is as follows:
%   f_R(r) = r/(sigma_x * sigma_y) exp(-ar^2) I_0(br^2)
%
% where sigma_x and sigma_y are the standard deviations of x and y, and
% a and b are defined:
%   a = (sigma_y^2 + sigma_x^2)/(2*sigma_x*sigma_y)^2
%   b = (sigma_y^2 - sigma_x^2)/(2*sigma_x*sigma_y)^2
%
% Nicholas O'Donoughue
% 23 Dec 2021

assert(all(rad>=0),'Ellipse radius CDF not defined for negative radii.');

% Eigendecomposition of covariance
lam = eig(cov); % each eigenvalue is the variance

if any(~isreal(lam))
%     warning('Complex eigenvalue encountered; the covariance matrix is likely slightly assymetric.');
    lam = eig(.5*(cov+cov'));
end

% Sort largest to smallest
lamSort = sort(lam,'descend');

%% Outlier Checking

% Check for degenerate condition; the ratio between the first two
% eigenvalues >> 1000
ratio = lamSort(1)/lamSort(2);

if isnan(ratio) || isinf(ratio) || ratio > 1e10
    % It's effectively a 1D problem.  We want to compute the Gaussian CDF
    % between -rad and rad
    
    % The largest eigenvalue is not finite; something screwy happened.
    % Let's set CDF to zero.
    if ~isfinite(lamSort(1))
        warning('Encountered a non-finite eigenvalue. Unable to compute ellipseCDF. Returning 0.');
        cdf = 0;
        return;
    end
    
    sigma = sqrt(lamSort(1));
    cdf = utils.normcdf(rad,0,sigma) - utils.normcdf(-rad,0,sigma);
    return
end

%% Pre-processing, Compute Parameters

% Extract variances and standard deviations
var_x = abs(lamSort(2));
var_y = abs(lamSort(1));
sigma_x = sqrt(var_x);
sigma_y = sqrt(var_y);

% Compute a and b
a = (var_y + var_x) / (4*var_x*var_y);
b = (var_y - var_x) / (4*var_x*var_y);

%% Build PDF numerical integration bounds
n_pts = 10001;
dr = rad/(n_pts-1);
r_vec = 0:dr:rad;

%% Evaluate PDF
pdf = (r_vec/(sigma_x*sigma_y)) .* exp((b-a).*r_vec.^2) .* besseli(0, b*r_vec.^2,1);
cdf = sum(pdf.*dr);