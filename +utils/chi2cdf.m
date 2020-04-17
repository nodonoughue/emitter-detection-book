function p = chi2cdf(x,v,tail)
% p = chi2cdf(x,v)
%
% Reduced functionality calculation of the CDF of a chi-squared random
% variable with v degrees of freedom.  To be used if the user does not
% have the Statistics & Machine Learning Toolbox installed.
%
% The optional third argument, tail, is a string with 'upper' or 'lower.'
% This specifies which end of the CDF is returned.  If unspecified, 'lower'
% will be used, corresponding to a conventional CDF (0 to x).
%
% When tail is upper, the built-in GAMMAINC function in MATLAB uses a more
% accurate version to compute the complement of the CDF than simply
% subtracting the CDF from 1.
%
% INPUTS:
%   x           The value at which to compute the CDF
%   v           Number of degrees of freedom
%   tail        [Optional] string specifying whether a normal CDF ('lower')
%               or complementary CDF ('upper') is desired.
%
% OUTPUTS:
%   p           Normal or complementary CDF for a Chi-Squared random
%               variable with v degrees of freedom, evaluated at x.
%
% Nicholas O'Donoughue
% 6 April 2020

if x < 0
    warning('Chi-Squared CDF requires a positive value of x; setting x=0.')
    x = 0;
end

if nargin < 2 || isempty(tail)
    tail = 'lower';
end

p = real(gammainc(x/2,v/2,tail));