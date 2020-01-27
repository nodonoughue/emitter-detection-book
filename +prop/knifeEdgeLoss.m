function L = knifeEdgeLoss(d1,d2,h)
% L = knifeEdgeLoss(d1,d2,h)
%
%
% Inputs:
%   d1      Distance from the transmitter to the obstruction [m]
%   d2      Distance from the obstruction to the receiver [m]
%   h       Vertical distance between the top of the obstruction and the
%           line of sight between the transmitter and receiver [m]
%
% Outputs:
%   L       Path loss [dB]
%
% Nicholas O'Donoughue
% 1 July 2019

% Equation (B.5)
nu = h.*sqrt(2)./(1+d1./d2);

% Initialize the output loss matrix
L = zeros(size(nu));

% First piece-wise component of (B.6)
L(nu<=0) = 0;

% Second piece-wise component of (B.6)
mask = nu > 0 & nu <= 2.4;
L(mask) = 6+9*nu(mask)-1.27*nu(mask).^2;

% Third piece-wise component of (B.6)
mask = nu > 2.4;
L(mask) = 13 + 20*log10(nu(mask));
