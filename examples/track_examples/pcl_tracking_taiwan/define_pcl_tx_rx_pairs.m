% Define PCL transmitters and receivers, define tx/rx pairs to process

%% Define Receivers
lat_rx = [22.660, 23.502, 24.384, 24.818];
lon_rx = [120.784,121.051, 121.216, 121.394];
alt_rx = geointerp(srtm, srtm_ref, lat_rx, lon_rx);

[rx_e, rx_n, rx_u] = geodetic2enu(lat_rx, lon_rx, alt_rx, ctr_lat, ctr_lon, ctr_up, wgs84Ellipsoid, angle_unit);

xrx = [rx_e(:), rx_n(:), rx_u(:)]'; % 3 x M

%% Transmitter Specs
tx_specs = struct('erp',10*log10(fm_high.power_kw(:))+30+3,...
    'freq',fm_high.freq_mhz(:)*1e6,...
    'bw',200e3,...
    't_int',.1);
rx_specs = struct('gain',10,...
    'az_beamwidth',30,...
    'el_beamwidth',10,...
    'bw',200e3,...
    'nf',5,...
    'loss',6);
rcs = 0;

%% Generate tx/rx assignment pairs

tx_ref_idx = [1,2,3,8,4,6,7,10,5,11,9,17,14,18,12,22];
rx_ref_idx = [ones(1,4), 2*ones(1,4), 3*ones(1,4), 4*ones(1,4)]; % 4 channels for each rx
ref_idx = [tx_ref_idx;rx_ref_idx];


%% Plot assignments
fig1=figure;
ax=worldmap(latlim, lonlim);
land = shaperead('landareas','UseGeoCoords',true);
hdl=geoshow(ax,land,'FaceColor',[.5, .7, .5],'EdgeColor','black');

hdl_tx=scatterm(fm_high.lat(tx_ref_idx),fm_high.lon(tx_ref_idx),30,[.7,.4,.4],'o','filled','DisplayName','Transmitters');
hdl_rx=scatterm(lat_rx,lon_rx,20,[.4,.4,.7],'s','filled','DisplayName','Receivers');


colors=get(groot,'DefaultAxesColorOrder');
for idx = 1:size(ref_idx,2)
    tx_idx = tx_ref_idx(idx);
    rx_idx = rx_ref_idx(idx);
    
    this_tx_lat = fm_high.lat(tx_idx);
    this_tx_lon = fm_high.lon(tx_idx);
    this_rx_lat = lat_rx(rx_idx);
    this_rx_lon = lon_rx(rx_idx);
    
    plotm([this_tx_lat,this_rx_lat],[this_tx_lon, this_rx_lon],'Color',colors(rx_idx,:));
end
utils.exportPlot(fig1,[prefix '_assignments']);

