function y = modulo2pi(x)
% Perform a 2*pi modulo operation, but with the result centered on zero, spanning
% from -pi to pi, rather than on the interval 0 to 2*pi.
%
% Input:
%   x       Input, any size
%
% Output:
%   y       Same size as x, modulo onto the interval -pi to pi
%
% 24 April 2023
% Nicholas O'Donoughue

% Shift the input so that zero is now pi
x_shift = x + pi;

% Perform a modulo operation; result is on the interval [0, 2*pi)
x_modulo = mod(x_shift, 2*pi);

% Undo the shift, so that a zero input is now a zero output.
% Result is on the interval [-pi, pi)
y = x_modulo - pi;
