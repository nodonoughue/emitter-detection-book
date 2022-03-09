% figs = book2_ex7_4()
%
% Executes Example 7.4 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle
%
% Nicholas O'Donoughue
% 11 January 2022

fprintf('Example 7.4...\n');

% Define sensor positions
x_aoa = [-5e3, 5e3;
         0, 0];
v_aoa = [0, 0;
         1000, 1000];
x_tgt = [ 3e3;
         25e3];

% Define sensor accuracy
num_dims = size(x_aoa,1);
num_sensors = size(x_aoa,2);
sigma_theta = 10;
sigma_psi = sigma_theta*pi/180;
C_df = sigma_psi^2*eye(num_sensors);
L = chol(C_df,'lower');

% Define pulse timing
pri = 1e-3;
T = 30; % observation period
num_pulses = floor(T/pri)+1;

% Define outputs
x_est = zeros(size(x_tgt,1),num_pulses);
cep = zeros(1,num_pulses);
x_init = [1; 1]*1e3;

% Outputs
zeta_full = zeros(num_dims,num_pulses);
x_aoa_full = zeros(num_dims,num_sensors,num_pulses);
ellipses = cell(1,num_pulses);
lobs = cell(1,num_pulses);

% Step through pulses
for idx=1:num_pulses
    % Update positions
    this_x_aoa = x_aoa + v_aoa * (idx-1) * pri;

    % Update function handles
    z_fun = @(x) triang.measurement(this_x_aoa, x);
    H_fun = @(x) triang.jacobian(this_x_aoa, x)';

    % Generate noisy measurements
    psi = z_fun(x_tgt);
    zeta = psi + L*randn(num_dims,1);


    if idx==1
        % Initialization
        this_x = triang.lsSoln(this_x_aoa, zeta, C_df, x_init);
        this_P = triang.computeCRLB(this_x_aoa, this_x, C_df);
    else
        % EKF Update
        [this_x, this_P] = tracker.ekfUpdate(prev_x, prev_P, zeta, C_df, z_fun, H_fun);
    end

    % Store the results and update the variables
    x_est(:, idx) = this_x;
    cep(idx) = utils.computeCEP50(this_P);

    % Store position and msmt
    x_aoa_full(:,:,idx) = this_x_aoa;
    zeta_full(:,idx) = zeta;
    ellipses{idx} = utils.drawErrorEllipse(x_tgt,this_P,101);

    zeta_plus = zeta + sigma_psi;
    zeta_minus = zeta - sigma_psi;

    lob = triang.drawLob(this_x_aoa, zeta, x_tgt,2);
    lob_plus = triang.drawLob(this_x_aoa, zeta_plus, x_tgt,2);
    lob_minus = triang.drawLob(this_x_aoa, zeta_minus, x_tgt,2);

    lob_fill = cat(2,lob_minus, fliplr(lob_plus), lob_minus(:,1,:));
    lobs{idx} = struct('lob',lob,'lob_fill',lob_fill);

    prev_x = this_x;
    prev_P = this_P;
end

%% Compute Errors
err = sqrt(sum(abs(x_est - x_tgt).^2,1));
time_vec = pri*(1:num_pulses);

%% Find Time to Achieved Desired Error
desired_cep = 100;
first_good_sample = find(cep < desired_cep,1,'first');
if isempty(first_good_sample)
    fprintf('More than %.2f s required to achieve %.2f m CEP50.\n', max(time_vec), desired_cep);
else
    fprintf('%.2f s required to achieve %.2f m CEP50.\n', time_vec(first_good_sample), desired_cep);
end


%% Make a Movie
utils.resetPlotSettings; % default settings don't work well with subplots
fig1=figure;
fig1.WindowState='maximized';

% Measurements
subplot(221);
hdl_msmt=plot(time_vec,zeta_full*180/pi,'DisplayName','Measurements');
hold on;
hdl_marker_msmt=scatter(time_vec(1),zeta_full(:,1)*180/pi,'ko','filled');
grid on;
legend;
xlabel('Time [s]');
ylabel('AOA [deg]');
legend('Location','NorthEast');
utils.excludeFromLegend(hdl_msmt(2:end));
utils.excludeFromLegend(hdl_marker_msmt);

% Error
subplot(222);
semilogy(time_vec,err,'DisplayName','Measured');
hold on;
plot(time_vec,cep,'DisplayName','Predicted (CEP_{50})');
hdl_marker_err=scatter(time_vec(1),err(1),'ko','filled');
hdl_marker_cep=scatter(time_vec(1),cep(1),'ko','filled');
xlabel('Time [s]');
ylabel('Error [m]');
legend('Location','NorthEast');
grid on;
utils.excludeFromLegend(hdl_marker_err);
utils.excludeFromLegend(hdl_marker_cep);

% XY
subplot(223);
num_tail = 50;
hdl_sensor=scatter(x_aoa_full(1,:,num_tail)/1e3,x_aoa_full(2,:,num_tail)/1e3,'^','filled','DisplayName','Sensors');
hold on;
alpha_vec = linspace(0,1,num_tail);

colors = get(0,'DefaultAxesColorOrder');
fill_color = colors(6,:);
this_lob = lobs{num_tail};
clear hdl_lob;
hdl_lob(num_sensors) = hdl_sensor; % initialize array with a graphics handle.
for idx=1:num_sensors
    hdl_lob(idx) = fill(squeeze(this_lob.lob_fill(1,:,idx)/1e3), squeeze(this_lob.lob_fill(2,:,idx))/1e3, fill_color,'FaceAlpha',.3,'EdgeColor','none','DisplayName','Measurements');    
end
utils.excludeFromLegend(hdl_lob(2:end));
hdl_est = scatter(x_est(1,1:num_tail)/1e3,...
                  x_est(2,1:num_tail)/1e3,'blue','filled',...
                    'AlphaData',alpha_vec,...
                    'MarkerFaceAlpha','flat',...
                    'DisplayName','Estimated Position');
this_ellipse = ellipses{num_tail};
hdl_ell = plot(this_ellipse(1,:)/1e3, this_ellipse(2,:)/1e3,'-.','DisplayName','Error Ellipse');
scatter(x_tgt(1)/1e3,x_tgt(2)/1e3,'s','filled','DisplayName','Target');
grid on;
legend('Location','NorthWest');
xlabel('X [km]');
ylabel('Y [km]');
xlim([-6,6]);
ylim([0,30]);

% XY Zoom
subplot(224);
alpha_vec = linspace(0,1,num_tail);

hdl_est_zoom = scatter(x_est(1,1:num_tail)/1e3,...
                  x_est(2,1:num_tail)/1e3,'blue','filled',...
                    'AlphaData',alpha_vec,...
                    'MarkerFaceAlpha','flat',...
                    'DisplayName','Estimated Position');
hold on;
hdl_ell_zoom = plot(this_ellipse(1,:)/1e3, this_ellipse(2,:)/1e3,'-.','DisplayName','Error Ellipse');
scatter(x_tgt(1)/1e3,x_tgt(2)/1e3,'s','filled','DisplayName','Target');
grid on;
legend('Location','NorthWest');
xlabel('X [km]');
ylabel('Y [km]');
offset = max(max(this_ellipse,[],2)-min(this_ellipse,[],2));
xlim(x_tgt(1)/1e3 + offset*[-1 1]/1e3);
ylim(x_tgt(2)/1e3 + offset*[-1 1]/1e3);


samples_per_frame = 25; % new frame every x samples
num_frames = 1+ceil((num_pulses-num_tail)/samples_per_frame);

for idx_frame = 1:num_frames
    % Find the current sample
    K = min(num_pulses, num_tail+samples_per_frame*(idx_frame-1));

    % Pull the data
    this_time = time_vec(K);
    this_x_aoa = x_aoa_full(:,:,K);
    this_zeta = zeta_full(:,K);
    this_err = err(K);
    this_cep = cep(K);
    this_x_est = x_est(:,K+(-num_tail+1:0));
    this_ellipse = ellipses{K};
    this_lob = lobs{K};

    % Update Sensor Positions
    hdl_sensor.XData = this_x_aoa(1,:)/1e3;
    hdl_sensor.YData = this_x_aoa(2,:)/1e3;
    
    % Update Measurement Marker
    for idx=1:num_sensors
        hdl_marker_msmt(idx).XData = this_time;
        hdl_marker_msmt(idx).YData = this_zeta(idx)*180/pi;
    end

    % Update Error Marker
    hdl_marker_err.XData = this_time;
    hdl_marker_err.YData = this_err;
    hdl_marker_cep.XData = this_time;
    hdl_marker_cep.YData = this_cep;

    % Update Measurements (x/y)
    for idx=1:num_sensors
        hdl_lob(idx).XData = this_lob.lob_fill(1,:,idx)/1e3;
        hdl_lob(idx).YData = this_lob.lob_fill(2,:,idx)/1e3;
    end

    % Update estimates
    hdl_est.XData = this_x_est(1,:)/1e3;
    hdl_est.YData = this_x_est(2,:)/1e3;
    hdl_est_zoom.XData = this_x_est(1,:)/1e3;
    hdl_est_zoom.YData = this_x_est(2,:)/1e3;

    % Update error ellipse
    hdl_ell.XData = this_ell(1,:)/1e3;
    hdl_ell.YData = this_ell(2,:)/1e3;
    hdl_ell_zoom.XData = this_ell(1,:)/1e3;
    hdl_ell_zoom.YData = this_ell(2,:)/1e3;

    % Redraw
    drawnow;
    pause(.001);

end