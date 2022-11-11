% figs=book2_ex8_2()
%
% Executes Example 8.2 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle
%
% Nicholas O'Donoughue
% 25 January 2022

fprintf('Example 8.2...\n');

%% Define ship and sensor positions
t_inc = 15;
t_max = 15*60;
t_vec = 0:t_inc:t_max; % seconds
num_time = numel(t_vec);

% Origin is 0,0.  Alt is 10 km. Velocity is +y at 100 m/s.
% Let's make it a function for simplicity
x_aoa_full = [0;0;10e3] + [0;100;0] * t_vec;

% Origin is 50 km, 50 km. Alt is 0 m.  Velocity is -x at 5 m/s.
ship_accel_power = .05;
a_tgt_full = cat(1,-ship_accel_power + 2*ship_accel_power*rand(2,num_time),zeros(1,num_time));
v_tgt_full = [-20;0;0] + cumsum(a_tgt_full*t_inc,2);
x_tgt_full = [50e3;50e3;0] + cumsum(v_tgt_full*t_inc,2);

%% Measurement Statistics
sigma_theta = .2; % deg
sigma_psi = sigma_theta*pi/180; % rad
num_msmt=2; % az/el
R = sigma_psi^2*eye(num_msmt); % measurement error covariance

%% Generate Measurements
do2daoa = num_msmt==2;
z = zeros(num_msmt,num_time);
for idx=1:num_time
    z(:,idx) = triang.measurement(x_aoa_full(:,idx),...
                                x_tgt_full(:,idx),do2daoa);
end
noise = sigma_psi * randn(num_msmt, num_time);
zeta = z + noise;

%% Set Up Tracker
% Track in x/y, even though a/c and tgt in x/y/z
sigma_a = .05;
num_dims = 3; % number of dimensions to use in state

[f_fun, q_fun, state_space] = tracker.makeKinematicModel('cv',num_dims,sigma_a^2);
num_states = state_space.num_states;
pos_idx = state_space.pos_idx;
vel_idx = state_space.vel_idx;
F = f_fun(t_inc); % generate state transition matrix
Q = q_fun(t_inc); % generate process noise covariance matrix

%% Initialize Track State
x_pred = zeros(num_states,1);
P_pred = zeros(num_states, num_states);

% Initialize position with AOA estimate from first three measurement
x_aoa_init = x_aoa_full(:,1:3);
z_init = reshape(z(:,1:3)',[],1); % collect the 3 az, then the 3 el msmts
C_init = kron(R,eye(3));

x_init_guess = [10e3;10e3;0];
[bnd, bnd_grad] = utils.constraints.fixedAlt(0,'flat');
x_init = triang.gdSolnFixed(x_aoa_init,z_init,C_init,x_init_guess,bnd,100,[],[],1e3,[],[],[]);
P_init = triang.computeCRLBfixed(x_aoa_init, x_init, C_init, bnd_grad, do2daoa);
% x_init = triang.gdSoln(x_aoa_init,z_init,C_init,x_init_guess);
% P_init = triang.computeCRLB(x_aoa_init,x_init,C_init);

x_pred(pos_idx) = x_init(1:num_dims);
P_pred(pos_idx, pos_idx) = P_init(1:num_dims,1:num_dims);

% Bound initial velocity uncertainty by assumed max velocity of 10 m/s
% (~20 knots)
max_vel = 30;
P_pred(vel_idx, vel_idx) = max_vel^2*eye(num_dims);
P_pred(vel_idx(end),:) = 0;
P_pred(:,vel_idx(end),:) = 0;

%% Step Through Time
fprintf('Iterating through EKF tracker time steps...');
x_ekf_est = zeros(num_dims,num_time);
x_ekf_pred = zeros(num_dims, num_time);
rmse_cov_est = zeros(1,num_time);
rmse_cov_pred = zeros(1, num_time);
num_ell_pts = 101; % number of points for ellipse drawing
x_ell_est = zeros(2,num_ell_pts,num_time);
lobs = cell(1,num_time);

for idx=1:num_time
    fprintf('.');

    % Grab Current Measurement
    this_zeta = zeta(:,idx);
    
    % Update msmt function
    this_x_aoa = x_aoa_full(1:num_dims,idx);
    this_x_tgt = x_tgt_full(1:num_dims,idx);

    [z_fun, h_fun] = tracker.makeMeasurementModel(this_x_aoa,[],[],[],[],[],state_space);

    % Update Position Estimate
    % Previous prediction stored in x_pred, P_pred
    % Updated estimate will be stored in x_est, P_est
    [x_est, P_est] = tracker.ekfUpdate(x_pred, P_pred, this_zeta, R, z_fun, h_fun);
    
    % Enforce known altitude
    x_est(pos_idx(end)) = 0;
    x_est(vel_idx(end)) = 0;
    P_est(pos_idx(end),:) = 0;
    P_est(:,pos_idx(end)) = 0;
    P_est(vel_idx(end),:) = 0;
    P_est(:,vel_idx(end)) = 0;
    
    % Predict state to the next time step
    [x_pred, P_pred] = tracker.kfPredict(x_est, P_est, Q, F);
    x_pred(pos_idx(end)) = 0;
    x_pred(vel_idx(end)) = 0;
    P_pred(pos_idx(end),:) = 0;
    P_pred(:,pos_idx(end)) = 0;
    P_pred(vel_idx(end),:) = 0;
    P_pred(:,vel_idx(end)) = 0;

    % Output the current prediction/estimation state
    x_ekf_est(:,idx) = x_est(pos_idx);
    x_ekf_pred(:,idx) = x_pred(pos_idx);

    rmse_cov_est(:,idx) = sqrt(trace(P_est(pos_idx,pos_idx)));
    rmse_cov_pred(:,idx)= sqrt(trace(P_pred(pos_idx,pos_idx)));

    % Draw an error ellipse
    x_ell_est(:,:,idx) = utils.drawErrorEllipse(x_est(pos_idx(1:2)), P_est(pos_idx(1:2), pos_idx(1:2)), num_ell_pts);

    % Measurement LOBs
    zeta_plus = this_zeta + sigma_psi;
    zeta_minus = this_zeta - sigma_psi;

    lob = triang.drawLob(this_x_aoa(1:2), this_zeta(1), this_x_tgt(1:2),5);
    lob_plus = triang.drawLob(this_x_aoa(1:2), zeta_plus(1), this_x_tgt(1:2),5);
    lob_minus = triang.drawLob(this_x_aoa(1:2), zeta_minus(1), this_x_tgt(1:2),5);

    lob_fill = cat(2,lob_minus, fliplr(lob_plus), lob_minus(:,1,:));
    lobs{idx} = struct('lob',lob,'lob_fill',lob_fill);
end

%% Compute Error
err_pred = x_ekf_pred(:,1:end-1) - x_tgt_full(1:num_dims,2:end);
err_est = x_ekf_est - x_tgt_full(1:num_dims,:);

rmse_pred = sqrt(sum(abs(err_pred).^2,1));
rmse_est = sqrt(sum(abs(err_est).^2,1));




%% Make a Movie
utils.resetPlotSettings; % default settings don't work well with subplots
fig1=figure;
fig1.WindowState='maximized';

% Measurements
subplot(221);
hdl_msmt=plot(t_vec,zeta*180/pi,'DisplayName','Measurements');
hold on;
hdl_marker_msmt=scatter(t_vec(1),zeta(:,1)*180/pi,'ko','filled');
grid on;
legend;
xlabel('Time [s]');
ylabel('AOA [deg]');
legend('Location','NorthEast');
utils.excludeFromLegend(hdl_msmt(2:end));
utils.excludeFromLegend(hdl_marker_msmt);

% Error
subplot(222);
semilogy(t_vec,rmse_est,'DisplayName','Measured');
hold on;
plot(t_vec,rmse_cov_est,'DisplayName','Predicted (CEP_{50})');
hdl_marker_err=scatter(t_vec(1),rmse_est(1),'ko','filled');
hdl_marker_cep=scatter(t_vec(1),rmse_cov_est(1),'ko','filled');
xlabel('Time [s]');
ylabel('Error [m]');
legend('Location','NorthEast');
grid on;
utils.excludeFromLegend(hdl_marker_err);
utils.excludeFromLegend(hdl_marker_cep);

% XY
subplot(223);
num_tail = 10;
hdl_sensor=scatter(x_aoa_full(1,num_tail)/1e3,x_aoa_full(2,num_tail)/1e3,'^','filled','DisplayName','Sensors');
hold on;

colors = get(0,'DefaultAxesColorOrder');
fill_color = colors(6,:);
this_lob = lobs{num_tail};
hdl_lob = fill(this_lob.lob_fill(1,:)/1e3, this_lob.lob_fill(2,:)/1e3, fill_color,'FaceAlpha',.3,'EdgeColor','none','DisplayName','Measurements');    

alpha_vec = linspace(0,1,num_tail);

hdl_est = scatter(x_ekf_est(1,1:num_tail)/1e3,...
                  x_ekf_est(2,1:num_tail)/1e3,'blue','filled',...
                    'AlphaData',alpha_vec,...
                    'MarkerFaceAlpha','flat',...
                    'DisplayName','Estimated Position');
this_ellipse = squeeze(x_ell_est(:,:,num_tail));
hdl_ell = plot(this_ellipse(1,:)/1e3, this_ellipse(2,:)/1e3,'-.','DisplayName','Error Ellipse');
hdl = plot(x_tgt_full(1,:)/1e3,x_tgt_full(2,:)/1e3);
utils.excludeFromLegend(hdl);
hdl_tgt = scatter(x_tgt_full(1,num_tail)/1e3,x_tgt_full(2,num_tail)/1e3,'s','filled','MarkerFaceColor',hdl.Color,'DisplayName','Target');
grid on;
legend('Location','SouthEast');
xlabel('X [km]');
ylabel('Y [km]');
xlim([-5,55]);
ylim([-5,105]);

% XY Zoom
subplot(224);
alpha_vec = linspace(0,1,num_tail);

hdl_est_zoom = scatter(x_ekf_est(1,1:num_tail)/1e3,...
                  x_ekf_est(2,1:num_tail)/1e3,'blue','filled',...
                    'AlphaData',alpha_vec,...
                    'MarkerFaceAlpha','flat',...
                    'DisplayName','Estimated Position');
hold on;
hdl_ell_zoom = plot(this_ellipse(1,:)/1e3, this_ellipse(2,:)/1e3,'-.','DisplayName','Error Ellipse');
hdl = plot(x_tgt_full(1,:)/1e3,x_tgt_full(2,:)/1e3);
utils.excludeFromLegend(hdl);
hdl_tgt_zoom = scatter(x_tgt_full(1,num_tail)/1e3,x_tgt_full(2,num_tail)/1e3,'s','filled','MarkerFaceColor',hdl.Color,'DisplayName','Target');
grid on;
legend('Location','SouthEast');
xlabel('X [km]');
ylabel('Y [km]');
offset = 10e3;
xlim(x_tgt_full(1,num_tail)/1e3 + offset*[-1,1]/1e3);
ylim(x_tgt_full(2,num_tail)/1e3 + offset*[-1,1]/1e3);


samples_per_frame = 1; % new frame every x samples
num_frames = 1+ceil((num_time-num_tail)/samples_per_frame);

for idx_frame = 1:num_frames
    % Find the current sample
    K = min(num_time, num_tail+samples_per_frame*(idx_frame-1));

    % Pull the data
    this_time = t_vec(K);
    this_x_aoa = x_aoa_full(:,K);
    this_x_tgt = x_tgt_full(:,K);
    this_zeta = zeta(:,K);
    this_err = rmse_est(K);
    this_cep = rmse_cov_est(K);
    this_x_est = x_ekf_est(:,K+(-num_tail+1:0));
    this_ellipse = squeeze(x_ell_est(:,:,K));
    this_lob = lobs{K};

    % Update Sensor/Target Positions
    hdl_sensor.XData = this_x_aoa(1,:)/1e3;
    hdl_sensor.YData = this_x_aoa(2,:)/1e3;
    hdl_tgt.XData = this_x_tgt(1,:)/1e3;
    hdl_tgt.YData = this_x_tgt(2,:)/1e3;
    hdl_tgt_zoom.XData = this_x_tgt(1,:)/1e3;
    hdl_tgt_zoom.YData = this_x_tgt(2,:)/1e3;

    % Update Measurement Marker
    for idx=1:num_msmt
        hdl_marker_msmt(idx).XData = this_time;
        hdl_marker_msmt(idx).YData = this_zeta(idx)*180/pi;
    end

    % Update Error Marker
    hdl_marker_err.XData = this_time;
    hdl_marker_err.YData = this_err;
    hdl_marker_cep.XData = this_time;
    hdl_marker_cep.YData = this_cep;

    % Update estimates
    hdl_est.XData = this_x_est(1,:)/1e3;
    hdl_est.YData = this_x_est(2,:)/1e3;
    hdl_est_zoom.XData = this_x_est(1,:)/1e3;
    hdl_est_zoom.YData = this_x_est(2,:)/1e3;

    % Update error ellipse
    hdl_ell.XData = this_ellipse(1,:)/1e3;
    hdl_ell.YData = this_ellipse(2,:)/1e3;
    hdl_ell_zoom.XData = this_ellipse(1,:)/1e3;
    hdl_ell_zoom.YData = this_ellipse(2,:)/1e3;
    
    % Measurement LOBs
    hdl_lob.XData = this_lob.lob_fill(1,:)/1e3;
    hdl_lob.YData = this_lob.lob_fill(2,:)/1e3;

    % Rescale zoomed plot
    subplot(224);
    xlim(this_x_tgt(1)/1e3 + 2*offset*[-1,1]/1e3);
    ylim(this_x_tgt(2)/1e3 + 2*offset*[-1,1]/1e3);

    % Redraw
    drawnow;
    pause(.01);

end