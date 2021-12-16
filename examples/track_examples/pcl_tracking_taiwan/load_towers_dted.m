% Read in towers and DTED
fm = readtable('fm_taiwan.xlsx');
fm = sortrows(fm, 'power_kw', 'descend');

latlim = [20 30];
lonlim = [115 125];

latlim_zoom = [21 26];
lonlim_zoom = [119 124];

addpath('C:\Users\nodonoug\Documents\GitLab\emitter-detection-book');
prefix = ['figures' filesep datestr(date(),'yyyy-mm-dd') filesep];

if ~exist(prefix,'dir')
    mkdir(prefix)
end

%% Sort Towers by Power
mask_low = fm.power_kw < 1;
mask_med = fm.power_kw < 10 & ~mask_low;
mask_high = ~mask_low & ~mask_med;

fm_high = fm(mask_high,:);

%% Pull in Terrain Data
if ~exist('srtm','var') || ~exist('srtm_ref','var')
    srtm_readerpath='C:\Users\nodonoug\OneDrive - RAND Corporation\Reference Info\DataSets\dted\SRTM';
    addpath(srtm_readerpath);
    level=1; % 1 = low, 3 = high
    % srtm_basedir='C:\Users\nodonoug\OneDrive - RAND Corporation\Reference Info\DataSets\dted\SRTM';
    srtm_basedir = 'D:\PRELUDE\PreludeData\DTED\SRTM';
    srtm_dir = [srtm_basedir, filesep, sprintf('Level%d',level)];
    [srtm,srtm_ref] = read_srtm_data(latlim_zoom,lonlim_zoom,srtm_dir,level);
end

%% Get Height for each transmitter
[fm_high.ht_m] = geointerp(srtm,srtm_ref,fm_high.lat, fm_high.lon);

%% Convert Transmitters to ECEF
ctr_lat = 23.5;
ctr_lon = 121.0;
ctr_up = 0;
angle_unit='deg';
[tx_e, tx_n, tx_u] = geodetic2enu(fm_high.lat, fm_high.lon, fm_high.ht_m, ctr_lat, ctr_lon, ctr_up,wgs84Ellipsoid,angle_unit);

xtx = [tx_e(:), tx_n(:), tx_u(:)]'; % 3 x N

