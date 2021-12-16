% Test the ability of a PCL to track an incoming aero target
%
% Assume up to 3G of lateral divert, up to 2G of vertical divert
%
% Otherwise, the target will follow a Maneuvering Re-entry Vehicle (MaRV)
% Kinematic Model
%
% Tracking is done in a 9-element state vector
%   x/y/z, vx/vy/vx, ax/ay/az
%
% Acceleration is broken into three components
%   aT - thrust
%   aA - aerodynamic
%   aG - gravity
%
% We assume that aG is known perfectly
%
% This script considers a target at a number of positions, with PCL
% transmitters and receivers locations defined statically.
%
% For each position, the PCL geolocation error is computed, and that error
% estimate, based on the CRLB, is used to compute the steady-state tracker
% error that can be expected.
%
% This is not the true tracker error at any point, which would depend on
% dynamics of target trajectory and noise.  Rather, this is a notion of how
% well we can track (in the steady state) a target at various locations,
% based on how well an update at that location will perform for an
% established track.

warning('off','all')

%% Load FM Towers and DTED data
load_towers_dted;

%% Define PCL TX/Rx Specs and Pairs for Processing
define_pcl_tx_rx_pairs;

%% Define Target Area
grid_ctr= [-200, 100, 9.144]*1e3; % 300 km due West, 30 kft alt
max_offset=300e3;
n_pts=101;
v_tgt=[.8*300, 0, 0]; % due East at .8 Mach

%% Test Flags
plot_all_links=false;

%% Setup Tracker Parameters
update_rate = 60; % seconds
heading = 135; % deg E of N
latency = 10; % processing latency
maxG=3; % assumed maneuverability
num_states = 2; % pos/vel/accel
num_dims = 3;

% Make State Matrices
F = tracker.makeTransitionMatrix(update_rate, num_states, num_dims);
H_posvel = cat(2,eye(2*num_dims), zeros(2*num_dims,num_dims*(num_states-2)));      % x/y/z/vx/vy/vz measurement
H_pos = cat(2,eye(num_dims),zeros(num_dims,num_dims*(num_states-1)));  % x/y/z measurement

target_type = '3d-divert';
Q = tracker.makeCAProcessNoise(maxG, num_states, heading, 0, update_rate, target_type);
% R = crlb

        
%% Run CRLB Analysis
% Rng Only PCL, No Doppler, No AOA
for do_doppler = [false, true]
    for angle_dims = 0:2
        this_prefix = [prefix 'track_analysis_vs_' target_type];
        if do_doppler
            this_prefix = [this_prefix '_dop'];
            suffix = 'Range-Doppler PCL';
        else
            this_prefix = [this_prefix '_nodop'];
            suffix = 'Range-only PCL';
        end
        if angle_dims == 0
            this_prefix = [this_prefix '_noaoa'];
            suffix = [suffix ', No AOA'];
        elseif angle_dims ==1
            this_prefix = [this_prefix '_1daoa'];
            suffix = [suffix, ', 1D AOA'];
        else
            this_prefix = [this_prefix '_2daoa'];
            suffix = [suffix, ', 2D AOA'];
        end
        
        [figs,crlb,~,~,~] = analyze_geometry(xtx,xrx,ref_idx,tx_specs,rx_specs,rcs,grid_ctr,max_offset,n_pts,v_tgt,do_doppler,angle_dims,plot_all_links);
        close(figs);
        
        num_test_cases = size(crlb,3);
        
        %% Track Analysis
        if do_doppler
            H = H_posvel;
        else
            H = H_pos;
        end
        
        % Outputs
        rmse_3d = zeros(size(latency));
        rmse_2d = zeros(size(latency));
        rmse_proj = zeros(size(latency));
        
        % Coordinates for Seeker FOV Projection
        x_tgt = [0, 0, 0]'; % Assume target at origin
        proj_az = 90; % deg E of N
        proj_el = 60; % deg above horiz
        x_obs = [sin(proj_az)*cos(proj_el); cos(proj_az)*cos(proj_el); sin(proj_el)];
        
        progress_tic = 100; % number of iterations per tic
        progress_newline = 40 * progress_tic; % number of tics per newline
        fprintf('Running tracker analysis for each target position...\n');
        for idx = 1:num_test_cases
            if mod(idx,progress_tic)==1, fprintf('.'); end
            if mod(idx,progress_newline)==0, fprintf('\n'); end
            
            % Grab the current measurement errors
            R = squeeze(crlb(:,:,idx));
            
            % Predict track error using horizontal divert model (steady-state,
            % level flight)
            % Rng-only PCL, No AOA, No Doppler
            [~, Pe] = tracker.steadystateError(F,H,Q,R);
            
            % Analyze
            Pe_pos = Pe(1:3,1:3);
            rmse_3d(idx) = sqrt(sum(abs(diag(Pe_pos))));
            rmse_2d(idx) = sqrt(sum(abs(diag(Pe_pos(1:2,1:2)))));
            P_proj = utils.projectError(Pe_pos, x_tgt, x_obs);
            rmse_proj(idx) = sqrt(sum(abs(diag(P_proj))));
        end
        fprintf('done.\n');
        
        %% Generate grid in ENU, convert to LLA
        grid_e = grid_ctr(1) + max_offset*linspace(-1,1,n_pts);
        grid_n = grid_ctr(2) + max_offset*linspace(-1,1,n_pts);
        grid_u = grid_ctr(3);
        [E,N] = meshgrid(grid_e,grid_n);
        [gridlat,gridlon,gridalt]  = enu2geodetic(E,N,grid_u,ctr_lat,ctr_lon,ctr_up,wgs84Ellipsoid,'deg');
        
        do_geo_plots = false;
        if do_geo_plots
            % Initialize Plot and Underlay
            fig1=figure;caxis autoc
            ax=worldmap(latlim, lonlim);
            land = shaperead('landareas','UseGeoCoords',true);
            geoshow(ax,land,'FaceColor',[.5, .7, .5],'EdgeColor','black');
            
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
            
            % Plot Projected RMSE
            z = reshape(rmse_proj,size(E));
            hdl=geoshow(ax,gridlat,gridlon,z/1e3,'DisplayType','texturemap','FaceAlpha',.6);
            hdl.ZData(:)=0;
            utils.excludeFromLegend(hdl)
            caxis(ax,[0 10]);
            cb=colorbar(ax);
            legend([hdl_tx,hdl_rx],'Location','NorthWest');
            title(ax,sprintf('Projected RMSE [km] for Steady State Track, %s', suffix));
            
            utils.exportPlot(fig1,[this_prefix '_proj_rmse']);
            delete(hdl)
            
            % Plot 3D RMSE
            z = reshape(rmse_3d,size(E));
            hdl=geoshow(ax,gridlat,gridlon,z/1e3,'DisplayType','texturemap','FaceAlpha',.6);
            hdl.ZData(:)=0;
            utils.excludeFromLegend(hdl)
            caxis(ax,[0 10]);
            cb=colorbar(ax);
            legend([hdl_tx,hdl_rx],'Location','NorthWest');
            title(ax,sprintf('3D RMSE [km] for Steady State Track, %s', suffix));
            
            utils.exportPlot(fig1,[this_prefix '_3d_rmse']);
            delete(hdl)
            
            % Plot 2D RMSE
            z = reshape(rmse_2d,size(E));
            hdl=geoshow(ax,gridlat,gridlon,z/1e3,'DisplayType','texturemap','FaceAlpha',.6);
            hdl.ZData(:)=0;
            utils.excludeFromLegend(hdl)
            caxis(ax,[0 10]);
            cb=colorbar(ax);
            legend([hdl_tx,hdl_rx],'Location','NorthWest');
            title(ax,sprintf('2D RMSE [km] for Steady State Track, %s', suffix));
            
            utils.exportPlot(fig1,[this_prefix '_2d_rmse']);
            delete(hdl)
        else
            % Plot 3D RMSE
            fig1=figure;
            z = reshape(rmse_3d,size(E));
            imagesc(grid_e/1e3,grid_n/1e3,z/1e3);set(gca,'ydir','normal');
            hold on;
            hdl_tx=scatter(xtx(1,:)/1e3, xtx(2,:)/1e3,30,[.7,.4,.4],'o','filled','DisplayName','Transmitters');
            hdl_rx=scatter(xrx(1,:)/1e3, xrx(2,:)/1e3,20,[.4,.4,.7],'s','filled','DisplayName','Receivers');
            
            colorbar;
            legend([hdl_tx,hdl_rx],'Location','NorthWest');
            title(sprintf('3D RMSE [km] for Steady State Track, %s', suffix));
            xlabel('E [km]');
            ylabel('N [km]');
            utils.exportPlot(fig1,[this_prefix '_3d_rmse']);
            
            % Plot 2D RMSE
            fig2=figure;
            z = reshape(rmse_2d,size(E));
            imagesc(grid_e/1e3,grid_n/1e3,z/1e3);set(gca,'ydir','normal');
            hold on;
            hdl_tx=scatter(xtx(1,:)/1e3, xtx(2,:)/1e3,30,[.7,.4,.4],'o','filled','DisplayName','Transmitters');
            hdl_rx=scatter(xrx(1,:)/1e3, xrx(2,:)/1e3,20,[.4,.4,.7],'s','filled','DisplayName','Receivers');
            
            colorbar;
            legend([hdl_tx,hdl_rx],'Location','NorthWest');
            title(sprintf('2D RMSE [km] for Steady State Track, %s',suffix));
            xlabel('E [km]');
            ylabel('N [km]');
            utils.exportPlot(fig1,[this_prefix '_2d_rmse']);
            
            % Plot Projected RMSE
            fig3=figure;
            z = reshape(rmse_proj,size(E));
            hdl=imagesc(grid_e/1e3,grid_n/1e3,z/1e3);set(gca,'ydir','normal');
            hold on;
            hdl_tx=scatter(xtx(1,:)/1e3, xtx(2,:)/1e3,30,[.7,.4,.4],'o','filled','DisplayName','Transmitters');
            hdl_rx=scatter(xrx(1,:)/1e3, xrx(2,:)/1e3,20,[.4,.4,.7],'s','filled','DisplayName','Receivers');
            
            colorbar;
            legend([hdl_tx,hdl_rx],'Location','NorthWest');
            title(sprintf('Projected RMSE [km] for Steady State Track, %s', suffix));
            xlabel('E [km]');
            ylabel('N [km]');
            utils.exportPlot(fig1,[this_prefix '_proj_rmse']);
            
        end
        close all
    end
end

warning('on')