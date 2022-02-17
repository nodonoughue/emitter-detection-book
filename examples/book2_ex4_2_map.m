% This modification of Example 4.2 plots the output on a map of the Earth
% using the Mapping Toolbox
%
% Nicholas O'Donoughue
% 7 Feb 2022

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
[tdoa_ref_vec, tdoa_test_vec] = utils.parseReferenceSensor(ref_tdoa, n_tdoa);
[fdoa_ref_vec, fdoa_test_vec] = utils.parseReferenceSensor(ref_fdoa, n_fdoa);
ref_idx =  [1:2*n_aoa,      2*n_aoa + tdoa_ref_vec,  2*n_aoa + n_tdoa + fdoa_ref_vec];
test_idx = [nan(1,2*n_aoa), 2*n_aoa + tdoa_test_vec, 2*n_aoa + n_tdoa + fdoa_test_vec];

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
% Define Map Extent
latlim = [min(x_sensor_lla(1,:)-7), max(x_sensor_lla(1,:)+7)];
lonlim = [min(x_sensor_lla(2,:)-7), max(x_sensor_lla(2,:)+7)];

% Initialize Map
figure;
worldmap(latlim,lonlim);
setm(gca, 'FFaceColor', [.85 .85 .99])
geoshow('landareas.shp','FaceColor',[.3 .5 .3],'EdgeColor','k');

% Add Sensor Positions
plot3m([1;1]*x_sensor_lla(1,:),[1;1]*x_sensor_lla(2,:),[zeros(1,4);x_sensor_lla(3,:)],'k--');
plot3m(x_sensor_lla(1,:), x_sensor_lla(2,:), x_sensor_lla(3,:),'kv','MarkerFaceColor','k');
view(0,80);

% Add TDOA Results
[grid_lat, grid_lon, ~] = utils.enu2lla(x_grid(1,:), x_grid(2,:), x_grid(3,:), ref_lat, ref_lon, ref_alt);
geoshow(reshape(grid_lat,size(XX)), reshape(grid_lon,size(XX)), rmse_crlb/1e3,'DisplayType','texturemap');
colorbar;
colormap(flipud(utils.viridis));
caxis([0 20]);
title('RMSE [km] from CRLB of 4-sensor Satellite TDOA Geometry')

% %% Plot Up Component of Error
% err_up = reshape(sqrt(crlb(3,3,:)),size(XX));
% 
% fig=figure;
% imagesc(xx_grid/1e3,yy_grid/1e3,err_up/1e3);
% hold on;
% colorbar;
% plot(x_sensor_enu(1,:)/1e3,x_sensor_enu(2,:)/1e3,'ko');
% caxis([5,10]);
% grid on;
% xlabel('E [km]');
% ylabel('N [km]');
% 
