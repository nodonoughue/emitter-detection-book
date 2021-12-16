function [figs,crlb,rmse_2d,rmse_3d,cep50] = analyze_geometry(xtx,xrx,ref_idx,tx_specs,rx_specs,rcs,grid_ctr,max_offset,n_pts,v_tgt,do_doppler,angle_dims,plot_all_links)

addpath('C:\Users\nodonoug\Documents\GitLab\emitter-detection-book');

if isscalar(n_pts)
    n_pts = [n_pts, n_pts, 1];
end
if isscalar(max_offset)
    max_offset = [max_offset, max_offset, 0];
end

%% Build Grid
x_vec = grid_ctr(1) + max_offset(1)*linspace(-1,1,n_pts(1));
y_vec = grid_ctr(2) + max_offset(2)*linspace(-1,1,n_pts(2));
z_vec = grid_ctr(3) + max_offset(3)*linspace(-1,1,n_pts(3));

[xx,yy,zz] = meshgrid(x_vec,y_vec,z_vec);

x_tgt = [xx(:),yy(:),zz(:)]'; % 3 x n_tgt

%% Compute Sensor Accuracy
do_pos_dependent_snr = true;
if do_pos_dependent_snr
    % Compute SNR
    snr_db = pcl.bistaticSNR(xtx,xrx,x_tgt,tx_specs.erp,tx_specs.bw*tx_specs.t_int,tx_specs.freq,rx_specs.gain,rcs,rx_specs.bw,rx_specs.nf,rx_specs.loss,ref_idx);
        % nTxRxPair x nTgt [dB] input snr
    snr_lin = 10.^(snr_db/10);


    % Estimate range/vel accuracy with equations 6.20 and 6.21 from
    % Malanowski
    ar = .68; % rng est coeff.
    av = .85; % vel est coeff.
    aa = .7; % my uneducated guess -- malanowski doesn't go into this, so I don't have a good coefficient for this
    rng_res = utils.constants.c./tx_specs.bw; 
    vel_res = utils.constants.c./(tx_specs.freq.*tx_specs.t_int); 
    az_res = rx_specs.az_beamwidth;
    el_res = rx_specs.el_beamwidth;
    
    % Handle non-scalar resolutions (e.g. transmitter dependent)
    if numel(rng_res) > 1 && numel(rng_res) ~= size(snr_lin,1)
        rng_res = rng_res(ref_idx(1,:),:);
    end
    if numel(vel_res) > 1 && numel(vel_res) ~= size(snr_lin,1)
        vel_res = vel_res(ref_idx(1,:),:);
    end
    if numel(az_res) > 1 && numel(az_res) ~= size(snr_lin,1)
        az_res = az_res(ref_idx(1,:),:);
    end
    if numel(el_res) > 1 && numel(el_res) ~= size(snr_lin,1)
        el_res = el_res(ref_idx(1,:),:);
    end
    
    rng_err = ar * rng_res ./ sqrt(snr_lin); % TODO: handle non-scalar bw, t_int 
    vel_err = av * vel_res ./ sqrt(snr_lin); 
    az_err = aa * az_res ./ sqrt(snr_lin);
    el_err = aa * el_res ./ sqrt(snr_lin);

    % Convert from time to bistatic range, freq to bistatic range-rate
    var_range = rng_err.^2;
    var_rrate = vel_err.^2;
    deg2rad = pi/180;
    var_az = (az_err*deg2rad).^2;
    var_el = (el_err*deg2rad).^2;

    % Build covariance matrices
    if do_doppler
        if angle_dims==2
            cov = cell2mat(reshape(arrayfun(@(idx) diag([var_range(:,idx);var_rrate(:,idx);var_az(:,idx);var_el(:,idx)]),1:size(var_range,2),'UniformOutput',false),1,1,size(var_range,2)));
        elseif angle_dims==1
            cov = cell2mat(reshape(arrayfun(@(idx) diag([var_range(:,idx);var_rrate(:,idx);var_az(:,idx)]),1:size(var_range,2),'UniformOutput',false),1,1,size(var_range,2)));
        else
            cov = cell2mat(reshape(arrayfun(@(idx) diag([var_range(:,idx);var_rrate(:,idx)]),1:size(var_range,2),'UniformOutput',false),1,1,size(var_range,2)));
        end
    else
        if angle_dims==2
            cov = cell2mat(reshape(arrayfun(@(idx) diag([var_range(:,idx);var_az(:,idx);var_el(:,idx)]),1:size(var_range,2),'UniformOutput',false),1,1,size(var_range,2)));
        elseif angle_dims==1
            cov = cell2mat(reshape(arrayfun(@(idx) diag([var_range(:,idx);var_az(:,idx)]),1:size(var_range,2),'UniformOutput',false),1,1,size(var_range,2)));
        else
            cov = cell2mat(reshape(arrayfun(@(idx) diag(var_range(:,idx)),1:size(var_range,2),'UniformOutput',false),1,1,size(var_range,2)));
        end
    end
else
    snr_db = 12.8;
    snr_lin = 10.^(snr_db/10);

    % Estimate range/vel accuracy with equations 6.20 and 6.21 from
    % Malanowski
    ar = .68; % rng est coeff.
    av = .85; % vel est coeff.
    aa = .7; % my uneducated guess -- malanowski doesn't go into this, so I don't have a good coefficient for this
    rng_res = utils.constants.c./tx_specs.bw; 
    vel_res = utils.constants.c./(tx_specs.freq.*tx_specs.t_int);
    az_res = rx_specs.az_beamwidth;
    el_res = rx_specs.el_beamwidth;

    rng_err = ar * rng_res ./ sqrt(snr_lin);
    vel_err = av * vel_res ./ sqrt(snr_lin); 
    az_err = aa * az_res ./ sqrt(snr_lin);
    el_err = aa * el_res ./ sqrt(snr_lin);

    % Convert from time to bistatic range, freq to bistatic range-rate
    var_range = rng_err.^2;
    var_rrate = vel_err.^2;
    deg2rad = pi/180;
    var_az = (az_err*deg2rad).^2;
    var_el = (el_err*deg2rad).^2;
    
    % Build covariance matrices
    cov_rng = var_range * eye(size(xtx,2)*size(xrx,2));
    cov_vel = var_rrate * eye(size(xtx,2)*size(xrx,2));
    cov_az = var_az * eye(size(xtx,2)*size(xrx,2));
    cov_el = var_el * eye(size(xtx,2)*size(xrx,2));
    if do_doppler
        if angle_dims == 2
            cov = blkdiag(cov_rng,cov_vel,cov_az,cov_el);
        elseif angle_dims == 1
            cov = blkdiag(cov_rng,cov_vel,cov_az);
        else
            cov = blkdiag(cov_rng,cov_vel);
        end
    else
        if angle_dims == 2
            cov = blkdiag(cov_rng,cov_az,cov_el);
        elseif angle_dims == 1
            cov = blkdiag(cov_rng,cov_az);
        else
            cov = cov_rng;
        end
    end
end

figs = [];
%% Plot Individual Tx-Rx Link SNR, Range Error, Range-Rate Error
if plot_all_links
    for idx_tx = 3%1:size(xtx,2)
        for idx_rx = 1%1:size(xrx,2)
            fig1=figure;
            idx_msmt = idx_tx + size(xtx,2)*(idx_rx-1);
            subplot(221); % Bistatic SNR
            imagesc(x_vec/1e3,y_vec/1e3,reshape(snr_db(idx_msmt,:),size(xx)));
            hold on;
            plot(xtx(1,idx_tx)/1e3,xtx(2,idx_tx)/1e3,'k^','DisplayName','Transmitter');
            plot(xrx(1,idx_rx)/1e3,xrx(2,idx_rx)/1e3,'ko','DisplayName','Receiver');
            xlabel('X [km]');
            ylabel('Y [km]');
            title(sprintf('Bistatic SNR [dB], BT=%d dB',fix(10*log10(tx_specs.bw*tx_specs.t_int))));
    %         legend('Location','NorthEast');
            axis equal;
            colorbar;
            caxis([0 30])
            colormap(utils.viridis)

            subplot(222); % Bistatic Range Error
            imagesc(x_vec/1e3,y_vec/1e3,reshape(rng_err(idx_msmt,:),size(xx)));
            hold on;
            plot(xtx(1,idx_tx)/1e3,xtx(2,idx_tx)/1e3,'k^','DisplayName','Transmitter');
            plot(xrx(1,idx_rx)/1e3,xrx(2,idx_rx)/1e3,'ko','DisplayName','Receiver');
            xlabel('X [km]');
            ylabel('Y [km]');
            title('Bistatic Range Error [m]');
    %         legend('Location','NorthEast');
            axis equal;
            colorbar;
            caxis([0 1000])
            colormap(flipud(utils.viridis))

            subplot(223); % Bistatic Range Rate Error
            imagesc(x_vec/1e3,y_vec/1e3,reshape(vel_err(idx_msmt,:),size(xx)));
            hold on;
            plot(xtx(1,idx_tx)/1e3,xtx(2,idx_tx)/1e3,'k^','DisplayName','Transmitter');
            plot(xrx(1,idx_rx)/1e3,xrx(2,idx_rx)/1e3,'ko','DisplayName','Receiver');
            xlabel('X [km]');
            ylabel('Y [km]');
            title('Bistatic Range Rate Error [m/s]');
    %         legend('Location','NorthEast');
            axis equal;
            colorbar;
            caxis([0 10])
            colormap(flipud(utils.viridis))

            subplot(224); % Rx AOA Error
            imagesc(x_vec/1e3,y_vec/1e3,reshape(az_err(idx_msmt,:),size(xx)));
            hold on;
            plot(xtx(1,idx_tx)/1e3,xtx(2,idx_tx)/1e3,'k^','DisplayName','Transmitter');
            plot(xrx(1,idx_rx)/1e3,xrx(2,idx_rx)/1e3,'ko','DisplayName','Receiver');
            xlabel('X [km]');
            ylabel('Y [km]');
            title('Bistatic Azimuth Error [deg]');
    %         legend('Location','NorthEast');
            axis equal;
            colorbar;
            caxis([0 10])
            colormap(flipud(utils.viridis))
        end
    end
end
%% Compute CRLB
if ~do_doppler
    v_tgt = [];
end

crlb = pcl.computeCRLB(xtx,xrx,x_tgt,[],[],v_tgt,cov,ref_idx,angle_dims);
    % 3 x 3 x n_tgt or 6 x 6 x n_tgt

crlb_pos = crlb(1:3,1:3,:);
if do_doppler
    crlb_vel = crlb(4:6,4:6,:);
end

%% Plot RMSE (x/y)
% rmse = reshape(arrayfun(@(x) sqrt(trace(crlb_pos(:,:,x))), 1:size(crlb_pos,3)),size(xx));
rmse_2d = reshape(arrayfun(@(x) sqrt(abs(trace(crlb_pos(1:2,1:2,x)))),1:size(crlb_pos,3)),size(xx));
rmse_3d = reshape(arrayfun(@(x) sqrt(abs(trace(crlb_pos(:,:,x)))),1:size(crlb_pos,3)),size(xx));
cep50 = reshape(utils.computeCEP50(crlb_pos),size(xx));

fig2=figure;
imagesc(x_vec/1e3,y_vec/1e3,rmse_2d);
set(gca,'ydir','normal');
hold on;
[C,h]=contour(xx/1e3,yy/1e3,rmse_2d,[10,100,1000,10000],'k');
clabel(C,h);
plot(xtx(1,:)/1e3,xtx(2,:)/1e3,'k^','DisplayName','Transmitter');
plot(xrx(1,:)/1e3,xrx(2,:)/1e3,'ko','DisplayName','Receiver');
xlabel('X [km]');
ylabel('Y [km]');
title('2D RMSE (x-y) [m] for PCL Geolocation');
legend('Location','NorthEast');
axis equal;
colorbar;
caxis([0 1e4])
colormap(flipud(utils.viridis))

fig3=figure;
imagesc(x_vec/1e3,y_vec/1e3,rmse_3d);
set(gca,'ydir','normal');
hold on;
[C,h]=contour(xx/1e3,yy/1e3,rmse_3d,[10,100,1000,10000],'k');
clabel(C,h);
plot(xtx(1,:)/1e3,xtx(2,:)/1e3,'k^','DisplayName','Transmitter');
plot(xrx(1,:)/1e3,xrx(2,:)/1e3,'ko','DisplayName','Receiver');
xlabel('X [km]');
ylabel('Y [km]');
title('3D RMSE [m] for PCL Geolocation');
legend('Location','NorthEast');
axis equal;
colorbar;
caxis([0 1e4])
colormap(flipud(utils.viridis))

fig4=figure;
imagesc(x_vec/1e3,y_vec/1e3,cep50);
set(gca,'ydir','normal');
hold on;
[C,h]=contour(xx/1e3,yy/1e3,cep50,[10,100,1000,10000],'k');
clabel(C,h);
plot(xtx(1,:)/1e3,xtx(2,:)/1e3,'k^','DisplayName','Transmitter');
plot(xrx(1,:)/1e3,xrx(2,:)/1e3,'ko','DisplayName','Receiver');
xlabel('X [km]');
ylabel('Y [km]');
title('CEP_{50} (dominant two dimensions) [m] for PCL Geolocation');
legend('Location','NorthEast');
axis equal;
colorbar;
caxis([0 1e4])
colormap(flipud(utils.viridis))

figs = [figs, fig2, fig3, fig4];