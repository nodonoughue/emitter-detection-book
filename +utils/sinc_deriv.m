function y = sinc_deriv(x)
% y = sinc_deriv(x)
%
% Returns the derivative of sinc(x), which is given
%       y= (x * cos(x) - sin(x)) / x^2
% for x ~= 0.  When x=0, y=0.  The input is in radians.
%
% NOTE: The MATLAB sinc function is defined sin(pi*x)/(pi*x).  Its usage
% will be different.  For example, if calling
%           y = sinc(x)
% then the corresponding derivative will be
%           z = sinc_deriv(pi*x);
%
% Input:
%
%   x       Radians input, may have arbitrary dimensions.
%
% Output:
%
%   y       Derivative of sinc(x)
%
% Nicholas O'Donoughue
% 1 July 2019


mask = x~=0;

y = zeros(size(x));
xx = x(mask);
y(mask) = (xx.*cos(xx) - sin(xx))./xx.^2;