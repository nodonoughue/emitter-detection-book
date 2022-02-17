function fig = book2_ex4_2()
% fig=book2_ex4_2()
%
% Executes Example 4.2 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle
%
% Nicholas O'Donoughue
% 22 July 2021

% Set up sensors
x_sensor_lla = [26 + 52/60, 27 + 52/60, 23 + 19/60, 24 + 19/60;
                -(68 + 47/60), -(72 + 36/60), -(69 + 47/60), -(74 + 36/60),;
                500e3 * ones(1,4)]; % deg N, deg E, m

% Define sensor velocity in ENU
v_abs = 7.61*1e3; % m/s, based on 500 km orbital height above Earth
heading_deg = 60; % deg N from E
v_sensor_enu = v_abs * [sind(heading_deg); cosd(heading_deg); 0] * ones(1,size(x_sensor_lla,2));

% Center of observation area
ref_lat = 26; % deg Lat (N)
ref_lon = -71; % deg Lon (E)
ref_alt = 0; % alt (m)

[e_sensor, n_sensor, u_sensor] = utils.lla2enu(x_sensor_lla(1,:), ...
                                               x_sensor_lla(2,:), ... 
                                               x_sensor_lla(3,:),...
                                               ref_lat, ref_lon, ref_alt);
x_sensor_enu = [e_sensor; n_sensor; u_sensor];

err_aoa_deg = 3;
err_tdoa_s = 1e-5;
err_fdoa_Hz = 100;
f_source_Hz = 1e9;

% Build grid of positions within 500km of source position (ENU origin)
xx_offset = linspace(-500e3,500e3,1001);

xx_grid = xx_offset;
yy_grid = xx_offset;
zz_grid = 0;

[XX,YY] = ndgrid(xx_grid,yy_grid);
x_grid = [XX(:), YY(:), zz_grid*ones(numel(XX),1)]';
v_source_enu = [0;0;0];

%% Sensor Selection
n_tdoa = size(x_sensor_lla,2);
n_aoa = n_tdoa;
n_fdoa = n_aoa;
ref_tdoa = n_tdoa;
ref_fdoa = n_fdoa;

% Manually do a reference/text index set
% 1:n_aoa are AOA sensors
% n_aoa + (1:n_tdoa) are TDOA sensors
% n_aoa + n_tdoa + (1:n_fdoa) are FDOA sensors
% [tdoa_ref_vec, tdoa_test_vec] = utils.parseReferenceSensor(ref_tdoa, n_tdoa);
% [fdoa_ref_vec, fdoa_test_vec] = utils.parseReferenceSensor(ref_fdoa, n_fdoa);
% ref_idx =  [1:2*n_aoa,      2*n_aoa + tdoa_ref_vec,  2*n_aoa + n_tdoa + fdoa_ref_vec];
% test_idx = [nan(1,2*n_aoa), 2*n_aoa + tdoa_test_vec, 2*n_aoa + n_tdoa + fdoa_test_vec];

%% Error Covariance Matrix
cov_psi = (err_aoa_deg*pi/180)^2*eye(2*n_aoa); % rad^2
cov_r = (err_tdoa_s*utils.constants.c)^2*eye(n_tdoa); % m^2
cov_rr = (err_fdoa_Hz*utils.constants.c/f_source_Hz)^2*eye(n_fdoa); % m^2/s^2

cov_x = blkdiag(cov_psi, cov_r, cov_rr);

% Compute the CRLB
crlb = hybrid.computeCRLB(x_sensor_enu, x_sensor_enu, x_sensor_enu, v_sensor_enu - v_source_enu, ...
                          x_grid, cov_x, ref_tdoa, ref_fdoa);

% Compute and display the RMSE
rmse_crlb = reshape(arrayfun(@(i) sqrt(trace(crlb(:,:,i))), 1:size(crlb,3)), size(XX));

%% Plotting
fig=figure;
imagesc(xx_grid/1e3,yy_grid/1e3,rmse_crlb/1e3);
hold on;
colorbar;
plot(x_sensor_enu(1,:)/1e3,x_sensor_enu(2,:)/1e3,'ko');
caxis([5,10]);
grid on;
xlabel('E [km]');
ylabel('N [km]');
