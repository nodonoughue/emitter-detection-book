function d_ang = ang_diff(ang1,ang2)
% Compute the angle difference, in radians, between two angles
%
% Avoids overflow issues when the angle bisecter crosses from +pi to -pi.
%
% Nicholas O'Donoughue
% 6 Aug 2020

d_ang = mod((ang2 - ang1) + pi,2*pi)-pi;