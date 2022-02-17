function figs = book2_ex8_2()
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
% 31 January 2022

fprintf('Example 8.2...\n');

%% Define ship and sensor positions
t_inc = 30;
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

fig1=figure;
hdl=plot(x_aoa_full(1,:), x_aoa_full(2,:));
utils.excludeFromLegend(hdl);
hold on;
plot(x_aoa_full(1,end), x_aoa_full(2,end),'-^','Color',hdl.Color,'DisplayName','AOA Sensor Trajectory');
hdl=plot(x_tgt_full(1,:), x_tgt_full(2,:));
utils.excludeFromLegend(hdl);
plot(x_tgt_full(1,end),x_tgt_full(2,end),'-<','Color',hdl.Color,'DisplayName','Target Trajectory');

grid on;
legend('Location','NorthEast');

%% Measurement Statistics
sigma_theta = 1; % deg
sigma_psi = sigma_theta*pi/180; % rad
num_msmt=2; % az/el
R = sigma_psi^2*eye(num_msmt); % measurement error covariance

% Add DF overlays
colors = get(0,'DefaultAxesColorOrder');
fill_color = colors(6,:);
for idx_overlay = [1,2,3,fix(num_time/2),num_time]
    this_ac_pos = x_aoa_full(1:2,idx_overlay);
    this_tgt_pos = x_tgt_full(1:2,idx_overlay);

    psi_true = triang.measurement(this_ac_pos,this_tgt_pos,false);
    psi_err_plus = psi_true + sigma_psi;
    psi_err_minus = psi_true - sigma_psi;

%     lob = triang.drawLob(this_ac_pos, psi_true, this_tgt_pos,2);
    lob_plus = triang.drawLob(this_ac_pos, psi_err_plus, this_tgt_pos,2);
    lob_minus = triang.drawLob(this_ac_pos, psi_err_minus, this_tgt_pos,2);

    lob_fill = cat(2,lob_minus, fliplr(lob_plus), lob_minus(:,1));

    hdl = fill(lob_fill(1,:), lob_fill(2,:), fill_color,'FaceAlpha',.3,'EdgeColor','none','DisplayName','DF Error');
    if idx_overlay > 1
        utils.excludeFromLegend(hdl);
    end
end

xlim([-5e3,55e3]);
ylim([-5e3,105e3]);

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
for idx=1:num_time
    fprintf('.');

    % Grab Current Measurement
    this_zeta = zeta(:,idx);
    
    % Update msmt function
    this_x_aoa = x_aoa_full(1:num_dims,idx);

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
end

fprintf('done.\n');
plot(x_init(1), x_init(2), '+','DisplayName','Initial Position Estimate');
plot(x_ekf_est(1,:),x_ekf_est(2,:),'-','DisplayName','EKF (est.)');
plot(x_ekf_pred(1,:),x_ekf_pred(2,:),'-','DisplayName','EKF (pred.)');
grid on;

% Plot some error ellipses
%for idx=1:10:num_time
%    this_ell = squeeze(x_ell_est(:,:,idx));
%    hdl = patch(this_ell(1,:), this_ell(2,:),hdl_est.Color,'FaceAlpha',.2,'DisplayName','1$\sigma$ Error (est.)');
%    if idx~=1
%        utils.excludeFromLegend(hdl);
%    end
%end

%% Zoomed Plot on Target
fig2=figure;
hdl=plot(x_tgt_full(1,:), x_tgt_full(2,:));
hold on;
utils.excludeFromLegend(hdl);
plot(x_tgt_full(1,end),x_tgt_full(2,end),'-<','Color',hdl.Color,'DisplayName','Target Trajectory');
plot(x_init(1), x_init(2), '+','DisplayName','Initial Position Estimate');
plot(x_ekf_est(1,:),x_ekf_est(2,:),'-','DisplayName','EKF (est.)');
plot(x_ekf_pred(1,:),x_ekf_pred(2,:),'-','DisplayName','EKF (pred.)');
grid on;
xlim([30e3,50e3]);
legend('Location','NorthWest');
utils.setPlotStyle(gca,{'widescreen','equal'});

%% Compute Error
err_pred = x_ekf_pred(:,1:end-1) - x_tgt_full(1:num_dims,2:end);
err_est = x_ekf_est - x_tgt_full(1:num_dims,:);

rmse_pred = sqrt(sum(abs(err_pred).^2,1));
rmse_est = sqrt(sum(abs(err_est).^2,1));

fig3=figure;
plot(t_vec,rmse_cov_est,'DisplayName','RMSE (est. cov.)');  
hold on;
plot(t_vec(2:end), rmse_cov_pred(1:end-1),'DisplayName','RMSE (pred. cov)');
set(gca,'ColorOrderIndex',1);
plot(t_vec,rmse_est,'--','DisplayName','RMSE (est. act.)');
plot(t_vec(2:end), rmse_pred, '--','DisplayName','RMSE (pred. act.)');
grid on;
xlabel('Time [sec]');
ylabel('Error [m]');
legend('Location','NorthEast')
set(gca,'yscale','log');
utils.setPlotStyle(gca,{'widescreen'});



%% Return Figure Handles
figs = [fig1, fig2, fig3];