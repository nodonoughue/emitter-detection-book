function [e, n, u] = correctENU(eg, ng, ug)
% [e, n, u] = correct(eg, ng, ug)
%
% Correct the input ENU coordinates to account for the Earth's rotation.
% The inputs eg (East) and ng (North) are assumed to refer to distances 
% along the Earth's surface, rather than in the ENU cartesian plane.  
% Similarly, the up coordinate is assumed to refer to altitude 
% perpendicular to the Earth at the point defined by e and n, rather than 
% perpendicular at the local origin.
%
% The output coordinates represent the same point, but references in the
% lcoal cartesian plane ENU.
%
% Assumes a spherical Earth.
%
% Inputs:
%   eg          East coordinate, along the Earth's surface
%   ng          North coordinate, along the Earth's surface
%   ug          Up coordinate, local to the destination point
%
% Outputs:
%   e           East coordinate in local ENU
%   n           North coordinate in local ENU
%   u           Up coordinate in local ENU
%
% Nicholas O'Donoughue
% 24 November 2021

%% Parse Inputs
bearing = atan2(eg,ng); % degrees east of north
R_g = sqrt(eg.^2 + ng.^2);

%% Correct for Earth curvature
[R_en, u] = utils.reckonSphereENU(R_g,ug);

%% Compute East and North coordinates
e = R_en .* sin(bearing);
n = R_en .* cos(bearing);