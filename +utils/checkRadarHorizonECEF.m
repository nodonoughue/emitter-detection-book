function mask = checkRadarHorizonECEF(x1, x2)
% Check the radar horizon between two points, in ECEF coordinates (meters)
%
% Assumes a smooth Earth with 4/3 effective radius for refraction of RF
% waves.
%
% Returns true for all visible points
%
% Nicholas O'Donoughue
% 9 Feb 2022

%% Compute Distance between all pairs of points
dist = utils.rng(x1,x2);

% Find local altitude of all points
[~,~,alt1] = utils.ecef2lla(x1(1,:), x1(2,:), x1(3,:));
[~,~,alt2] = utils.ecef2lla(x2(1,:), x2(2,:), x2(3,:));

%% Compute radar horizon
horizon_range = prop.radarHorizon(alt1(:), alt2(:)'); % size(x1,2) x size(x2,2)

% Compare to straight line distance
mask = dist < horizon_range;