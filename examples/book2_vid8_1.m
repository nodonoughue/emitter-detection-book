% figs=book2_ex8_1()
%
% Executes Example 8.1 from Practical Geolocation for Electronic Warfare
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

fprintf('Example 8.1...\n');

% Define sensor positions
x_tdoa = [5e3,   0,  0, -5e3;
            0, 5e3,  0,    0;
           30,  60, 30,   60];
v_tdoa = 0; % all sensors stationary
[num_dims, n_tdoa] = size(x_tdoa);

%% Define target trajectory
x_tgt_init = [-50e3;
         100e3;
         20e3*unitsratio('m','ft')];
vel = 200;

t_e_leg = 3*60; % turn at 3 min
turn_rad = 50e3;
t_turn = pi/2*turn_rad/vel;
t_s_leg = 6*60;

t_inc = 10; % 10 seconds between track updates

% Due East
x_e_leg = x_tgt_init + [vel;0;0]*(0:t_inc:t_e_leg);

% Turn to South
angle_turn = pi/2 * (t_inc:t_inc:t_turn)/t_turn;
x_turn = x_e_leg(:,end) + turn_rad * [sin(angle_turn); cos(angle_turn)-1; zeros(size(angle_turn))];

% Due South
x_s_leg = x_turn(:,end) + [0;-vel;0]*(t_inc:t_inc:t_s_leg);

% Combine legs
x_tgt_full = cat(2,x_e_leg, x_turn, x_s_leg);
t_vec = 0:t_inc:(t_e_leg+t_turn+t_s_leg);
num_time = numel(t_vec);

%% Measurement Statistics
ref_idx = 1;
sigma_toa = 10e-9;
C_toa = sigma_toa^2 * eye(n_tdoa);
C_roa = utils.constants.c^2*C_toa;
R = utils.resampleCovMtx(C_roa, ref_idx);
L = chol(R,'lower');
num_msmt = size(R,1);

%% Generate Measurements
z = tdoa.measurement(x_tdoa, x_tgt_full, ref_idx);
noise = L * randn(num_msmt, num_time);
zeta = z + noise;

% crlb = tdoa.computeCRLB(x_tdoa, x_tgt_full, C_roa, ref_idx, false, true);
% rmse_crlb = sqrt(arrayfun(@(i) trace(crlb(:,:,i)), 1:size(crlb,3)));

%% Set Up Tracker
sigma_a = 10; % Try .1 and 10, to compare results

[f_fun, q_fun, state_space] = tracker.makeKinematicModel('cv',num_dims,sigma_a^2);
num_states = state_space.num_states;
pos_idx = state_space.pos_idx;
vel_idx = state_space.vel_idx;
F = f_fun(t_inc); % generate state transition matrix
Q = q_fun(t_inc); % generate process noise covariance matrix

[z_fun, h_fun] = tracker.makeMeasurementModel([],x_tdoa,[],v_tdoa,ref_idx,[],state_space);
 % msmt function and linearized msmt function

%% Initialize Track State
x_pred = zeros(num_states,1);
P_pred = zeros(num_states, num_states);

% Initialize position with TDOA estimate from first measurement
x_init = [0;50e3;5e3];
epsilon = 100;
x_pred(pos_idx) = tdoa.gdSoln(x_tdoa, zeta(:,1), C_roa, x_init,[],[],epsilon,[],[],[],ref_idx);
P_pred(pos_idx,pos_idx) = 10*tdoa.computeCRLB(x_tdoa,x_pred(pos_idx),C_roa,ref_idx,false,true);


% Bound initial velocity uncertainty by assumed max velocity of 340 m/s
% (Mach 1 at sea level)
max_vel = 340;
P_pred(vel_idx, vel_idx) = 10*max_vel^2*eye(num_dims);

%% Step Through Time
fprintf('Iterating through EKF tracker time steps...');
x_ekf_est = zeros(num_dims,num_time);
x_ekf_pred = zeros(num_dims, num_time);
rmse_cov_est = zeros(1,num_time);
rmse_cov_pred = zeros(1, num_time);
num_ell_pts = 101; % number of points for ellipse drawing
x_ell_est = zeros(2,num_ell_pts,num_time);
for idx=1:num_time
    fprintf('.');

    % Grab Current Measurement
    this_zeta = zeta(:,idx);

    % Update Position Estimate
    % Previous prediction stored in x_pred, P_pred
    % Updated estimate will be stored in x_est, P_est
    [x_est, P_est] = tracker.ekfUpdate(x_pred, P_pred, this_zeta, R, z_fun, h_fun);
    
    % Predict state to the next time step
    [x_pred, P_pred] = tracker.kfPredict(x_est, P_est, Q, F);

    % Output the current prediction/estimation state
    x_ekf_est(:,idx) = x_est(pos_idx);
    x_ekf_pred(:,idx) = x_pred(pos_idx);

    rmse_cov_est(:,idx) = sqrt(trace(P_est(pos_idx,pos_idx)));
    rmse_cov_pred(:,idx)= sqrt(trace(P_pred(pos_idx,pos_idx)));

    % Draw an error ellipse
    x_ell_est(:,:,idx) = utils.drawErrorEllipse(x_est(pos_idx(1:2)), P_est(pos_idx(1:2), pos_idx(1:2)), num_ell_pts);
end

fprintf('done.\n');

%% Compute Error
err_pred = x_ekf_pred(:,1:end-1) - x_tgt_full(:,2:end);
err_est = x_ekf_est - x_tgt_full;

rmse_pred = sqrt(sum(abs(err_pred).^2,1));
rmse_est = sqrt(sum(abs(err_est).^2,1));


%% Make a Movie
utils.resetPlotSettings; % default settings don't work well with subplots
fig1=figure;
fig1.WindowState='maximized';

% Measurements
subplot(221);
hdl_msmt=plot(t_vec,zeta,'DisplayName','Measurements');
hold on;
hdl_marker_msmt=scatter(t_vec(1),zeta(:,1),'ko','filled');
grid on;
legend;
xlabel('Time [s]');
ylabel('RDOA [m]');
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
scatter(x_tdoa(1,:)/1e3,x_tdoa(2,:)/1e3,'^','filled','DisplayName','Sensors');
hold on;
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
legend('Location','SouthWest');
xlabel('X [km]');
ylabel('Y [km]');

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
legend('Location','SouthWest');
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
    this_x_tgt = x_tgt_full(:,K);
    this_zeta = zeta(:,K);
    this_err = rmse_est(K);
    this_cep = rmse_cov_est(K);
    this_x_est = x_ekf_est(:,K+(-num_tail+1:0));
    this_ellipse = squeeze(x_ell_est(:,:,K));

    % Update Sensor Positions
    hdl_tgt.XData = this_x_tgt(1,:)/1e3;
    hdl_tgt.YData = this_x_tgt(2,:)/1e3;
    hdl_tgt_zoom.XData = this_x_tgt(1,:)/1e3;
    hdl_tgt_zoom.YData = this_x_tgt(2,:)/1e3;

    % Update Measurement Marker
    for idx=1:num_msmt
        hdl_marker_msmt(idx).XData = this_time;
        hdl_marker_msmt(idx).YData = this_zeta(idx);
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
    
    % Rescale zoomed plot
    subplot(224);
    xlim(this_x_tgt(1)/1e3 + 2*offset*[-1,1]/1e3);
    ylim(this_x_tgt(2)/1e3 + 2*offset*[-1,1]/1e3);

    % Redraw
    drawnow;
    pause(.001);

end