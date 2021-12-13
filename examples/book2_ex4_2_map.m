%% Define Map Extent
latlim = [min(x_sensor_lla(1,:)-10), max(x_sensor_lla(1,:)+10)];
lonlim = [min(x_sensor_lla(2,:)-10), max(x_sensor_lla(2,:)+10)];

%% Initialize Map
figure;
worldmap(latlim,lonlim);
setm(gca, 'FFaceColor', [.6 .6 .85])
geoshow('landareas.shp','FaceColor',[.3 .5 .3],'EdgeColor','k');

%% Add Sensor Positions
plot3m([1;1]*x_sensor_lla(1,:),[1;1]*x_sensor_lla(2,:),[zeros(1,4);x_sensor_lla(3,:)],'k--.');
view(0,30);

%% Add TDOA Results
[grid_lat, grid_lon, grid_alt] = utils.enu2lla(x_grid(1,:), x_grid(2,:), x_grid(3,:), ref_lat, ref_lon, ref_alt);
geoshow(reshape(grid_lat,size(XX)), reshape(grid_lon,size(XX)), rmse_crlb/1e3,'DisplayType','texturemap');
colorbar;
colormap(flipud(utils.viridis));
caxis([0 20]);

%% Plot Up Component of Error
err_up = reshape(sqrt(crlb(3,3,:)),size(XX));

fig=figure;
imagesc(xx_grid/1e3,yy_grid/1e3,e/1e3);
hold on;
colorbar;
plot(x_sensor_enu(1,:)/1e3,x_sensor_enu(2,:)/1e3,'ko');
caxis([5,10]);
grid on;
xlabel('E [km]');
ylabel('N [km]');

