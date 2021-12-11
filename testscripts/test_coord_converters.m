% Script to test forward/back conversion of coordinates
%
% Nicholas O'Donoughue
% 18 July 2021

%% ENU to AER

num_test_pts = 1000;

max_val = 1000e3; % 100 km max offset from local origin

% Generate a random position in ENU
pos_enu = max_val * (2*rand(3,num_test_pts)-1);
    
% Convert to AER
[az, el, range] = utils.enu2aer(pos_enu(1,:), pos_enu(2,:), pos_enu(3,:));

% Convert back to ENU
[e_out, n_out, u_out] = utils.aer2enu(az, el, range);

% Compute errors
pos_err = sqrt((e_out-pos_enu(1,:)).^2 + (n_out-pos_enu(2,:)).^2 ...
               + (u_out - pos_enu(3,:)).^2);

max_err = max(pos_err);
med_err = median(pos_err);
avg_err = mean(pos_err);

fprintf('Testing ENU-AER-ENU conversion for self-consistency...\n');
fprintf('%d random points tested with ENU coordinates <= %d km.\n',num_test_pts, max_val/1e3);
fprintf('Max Error: %.2g km, Median Error: %.2g km, Mean Error: %.2g km.\n',max_err/1e3, med_err/1e3, avg_err/1e3);

%% ENU to ECEF
% Use the same set of random ENU positions; but this time also assign a
% random lat/lon between +/- 60 lat.

lat0 = 60 * (2 * rand(1,num_test_pts) -1);
lon0 = 180 * (2 * rand(1,num_test_pts) - 1);

[x, y, z] = utils.enu2ecef(pos_enu(1,:), pos_enu(2,:), pos_enu(3,:), ...
    lat0, lon0, 0.0);

[e_out, n_out, u_out] = utils.ecef2enu(x, y, z, lat0, lon0, 0.0);

pos_err = sqrt((e_out-pos_enu(1,:)).^2 + (n_out-pos_enu(2,:)).^2 ...
               + (u_out - pos_enu(3,:)).^2);

max_err = max(pos_err);
med_err = median(pos_err);
avg_err = mean(pos_err);

fprintf('Testing ENU-ECEF-ENU conversion for self-consistency...\n');
fprintf('%d random points tested with ENU coordinates <= %d km.\n',num_test_pts, max_val/1e3);
fprintf('Max Error: %.2g km, Median Error: %.2g km, Mean Error: %.2g km.\n',max_err/1e3, med_err/1e3, avg_err/1e3);

%% ECEF to Geodetic

max_val = 8000e3;

pos_ecef = max_val * (2 * rand(3,num_test_pts) -1 );

[lat, lon, alt] = utils.ecef2lla(pos_ecef(1,:), pos_ecef(2,:), pos_ecef(3,:));

[x_out, y_out, z_out] = utils.lla2ecef(lat, lon, alt);

pos_err = sqrt((x_out - pos_ecef(1,:)).^2 + (y_out - pos_ecef(2,:)).^2 ...
               +(z_out - pos_ecef(3,:)).^2);

max_err = max(pos_err);
med_err = median(pos_err);
avg_err = mean(pos_err);

fprintf('Testing ECEF-LLA-ECEF conversion for self-consistency...\n');
fprintf('%d random points tested with ENU coordinates <= %d km.\n',num_test_pts, max_val/1e3);
fprintf('Max Error: %.2g km, Median Error: %.2g km, Mean Error: %.2g km.\n',max_err/1e3, med_err/1e3, avg_err/1e3);

%% ECEF to AER

% Define a lat/lon ref point for AER
lat0 = 60 * (2 * rand(1,num_test_pts) -1);
lon0 = 180 * (2 * rand(1,num_test_pts) - 1);

[az, el, range] = utils.ecef2aer(pos_ecef(1,:), pos_ecef(2,:), pos_ecef(3,:), ...
    lat0, lon0, 0.0);

[x_out, y_out, z_out] = utils.aer2ecef(az, el, range, lat0, lon0, 0.0);

pos_err = sqrt((x_out - pos_ecef(1,:)).^2 + (y_out - pos_ecef(2,:)).^2 ...
               +(z_out - pos_ecef(3,:)).^2);

max_err = max(pos_err);
med_err = median(pos_err);
avg_err = mean(pos_err);

fprintf('Testing ECEF-AER-ECEF conversion for self-consistency...\n');
fprintf('%d random points tested with ENU coordinates <= %d km.\n',num_test_pts, max_val/1e3);
fprintf('Max Error: %.2g km, Median Error: %.2g km, Mean Error: %.2g km.\n',max_err/1e3, med_err/1e3, avg_err/1e3);
