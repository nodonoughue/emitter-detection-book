function p = normcdf(x,mu,sigma)
% p = normcdf(x,mu,sigma)
%
% Compute the CDF for variable x, assuming a Gaussian
% distribution with mean mu and variance sigma.
%
% All non-singleton dimensions of inputs must be the same.
%
% INPUTS:
%   x       n-dimensional array of points at which to calculate CDF
%   mu      n-dimensional array of mean values for CDF
%   sigma   n-dimensional array of standard deviations for CDF
%
% OUTPUTS:
%   p       n-dimensional array of CDF outputs
%
% Nicholas O'Donoughue
% 14 December 2021

if nargin < 2 || isempty(mu)
    mu = 0;
end

if nargin < 3 || isempty(sigma)
    sigma = 1;
end

p = .5 * erfc(-(x-mu)./(sqrt(2)*sigma));