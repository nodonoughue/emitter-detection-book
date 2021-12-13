function [R_en, u] = reckonSphereENU(R_g, alt)
% [R_en, u] = reckonSphereENU(R_g, alt)
%
% Given some distance along the Earth's surface r_g, and and altitude above
% the Earth's surface alt, compute the in-plane distance parallel to the
% Earth's surface (local East-North) and the out-of-plane distance (up).
%
% Inputs:
%   R_g         Ground range, along the Earth's surface
%   alt         Altitude, above the Earth's surface
%
% Outputs:
%   R_en        Distance in the East-North plane
%   u           Distance perpendicular to the East-North plane
%
% Nicholas O'Donoughue
% 24 November 2021


%% Parse Inputs
Re = utils.constants.Re_true; % for simplicity

%% Find central angle - phi
phi = R_g / Re;

%% Radius and EN distance to the crossing of the radial and the EN plane
radial = Re ./ cos(phi);
en_dist = Re .* tan(phi);

%% Find the radial to get to the desired point
y = alt + Re - radial;

%% Compute the excess in-plane and the up
u = y.*cos(phi);
en_delta = y.*sin(phi);

%% Compute the distance in the 'east-north' plane
R_en = en_dist + en_delta;
