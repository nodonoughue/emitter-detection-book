function x = chi2inv(p,v)
% x = chi2inv(p,v)
%
% Utilize elementary functions to compute the inverse of the
% Chi-Squared CDF, for users that do not have the Statistics &
% Machine Learning Toolbox installed.
%
% This is achieved via the function gammaincinv, which computes
% the inverse of the incomplete Gamma function.
%
%
% Nicholas O'Donoughue
% 7 April 2020

x = gammaincinv(p,v/2)*2;
