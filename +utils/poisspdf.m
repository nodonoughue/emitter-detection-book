function y = poisspdf(x,lambda)
%POISSPDF Poisson probability density function.
%   Y = POISSPDF(X,LAMBDA) returns the Poisson probability density 
%   function with parameter LAMBDA at the values in X.
%
%   The size of Y is the common size of X and LAMBDA. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Note that the density function is zero unless X is an integer.
%
%   See also POISSCDF, POISSFIT, POISSINV, POISSRND, POISSTAT, PDF.
%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.22.
%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 2.13.2.5 $  $Date: 2004/12/06 16:38:01 $
if nargin <  2, 
    error('stats:poisspdf:TooFewInputs','Requires two input arguments.'); 
end
% [errorcode x lambda] = distchck(2,x,lambda);
% if errorcode > 0
%     error('stats:poisspdf:InputSizeMismatch',...
%           'Requires non-scalar arguments to match in size.');
% end
if ~isscalar(x)
    sz = size(x);
else
    sz = size(lambda);
end

if isa(x,'single') || isa(lambda,'single')
   y = zeros(sz,'single');
else
   y = zeros(sz);
end
if ~isfloat(x)
   x = double(x);
end
if (length(y) == 0), return; end
y(lambda < 0) = NaN;
k = (x >= 0 & x == round(x) & lambda > 0);
% using exp(gammaln) instead of gamma can avoid possible overflow.
if (any(k(:)))
   if ~isscalar(lambda)
       lambda_local = lambda(k);
   else
       lambda_local = lambda;
   end
   
   if ~isscalar(x)
       x_local = x(k);
   else
       x_local = x;
   end
   
   y(k) = exp(-lambda_local + x_local .* log(lambda_local) - gammaln(x_local + 1));
end
y(x==0 & lambda==0) = 1;