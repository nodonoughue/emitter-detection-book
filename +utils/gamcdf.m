function cdf = gamcdf (x, a, b)
% GAMCDF  CDF of the Gamma distribution
%  CDF = gamcdf(X, A, B) computes the cumulative distribution
%  function (CDF) at X of the Gamma distribution with parameters
%  A and B (i.e. mean of the distribution is A*B and variance
%  is A*B^2).

% Adapted for Matlab (R) from GNU Octave 3.0.1
% Original file: statistics/distributions/gamcdf.m
% Original author: TT <Teresa.Twaroch@ci.tuwien.ac.at>

% Copyright (C) 1995, 1996, 1997, 2005, 2006, 2007 Kurt Hornik
% Copyright (C) 2008-2009 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if (nargin ~= 3)
    error ('gamcdf: you must give three arguments');
end

sz = size (x);
cdf = zeros (sz);

k = find (~(a > 0) | ~(b > 0) | isnan (x));
if (any (k))
    cdf (k) = NaN;
end

k = find ((x > 0) & (a > 0) & (b > 0));
if (any (k))
    if ~isscalar(a)
        a = a(k);
    end
    
    if ~isscalar(b)
        b = b(k);
    end
    
    if ~isscalar(x)
        x = x(k);
    end
    
    cdf (k) = gammainc (x ./ b, a);
end

end

