function x = norminv(p,mu,sigma)
% x = norminv(p,mu,sigma)
%
% Compute the inverse CDF for probability p, assuming a Gaussian
% distribution with mean mu and variance sigma.
%
% All non-singleton dimensions of inputs must be the same.
%
% INPUTS:
%   p       n-dimensional array of CDF probabilities
%   mu      n-dimensional array of mean values for CDF
%   sigma   n-dimensional array of standard deviations for CDF
%
% OUTPUTS:
%   x       n-dimensional array of points at which to calculate CDF
%
% Nicholas O'Donoughue
% 14 December 2021

x = mu -sqrt(2) * erfcinv(2*p) .* sigma;