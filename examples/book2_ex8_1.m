function figs = book2_ex8_1()
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

% Plot Geometry
fig1 = figure;
plot(x_tdoa(1,:), x_tdoa(2,:), 'o', 'DisplayName','Sensors');
hold on;
hdl=plot(x_tgt_full(1,end), x_tgt_full(2,end), 'v-','DisplayName','Aircraft');
hdl_traj=plot(x_tgt_full(1,:),x_tgt_full(2,:),'Color',hdl.Color);
utils.excludeFromLegend(hdl_traj);
legend('Location','SouthWest');
grid on;

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
sigma_a = 1;

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
plot(x_init(1), x_init(2), '+','DisplayName','Initial Position Estimate');
hdl_est=plot(x_ekf_est(1,:),x_ekf_est(2,:),'--','DisplayName','EKF (est.)');
plot(x_ekf_pred(1,:),x_ekf_pred(2,:),'--','DisplayName','EKF (pred.)');
grid on;

% Plot some error ellipses
%for idx=1:10:num_time
%    this_ell = squeeze(x_ell_est(:,:,idx));
%    hdl = patch(this_ell(1,:), this_ell(2,:),hdl_est.Color,'FaceAlpha',.2,'DisplayName','1$\sigma$ Error (est.)');
%    if idx~=1
%        utils.excludeFromLegend(hdl);
%    end
%end
utils.setPlotStyle(gca,{'equal'});

%% Compute Error
err_pred = x_ekf_pred(:,1:end-1) - x_tgt_full(:,2:end);
err_est = x_ekf_est - x_tgt_full;

rmse_pred = sqrt(sum(abs(err_pred).^2,1));
rmse_est = sqrt(sum(abs(err_est).^2,1));

fig2=figure;
plot(t_vec,rmse_cov_est,'DisplayName','RMSE (est. cov.)');  
hold on;
plot(t_vec(2:end), rmse_cov_pred(1:end-1),'DisplayName','RMSE (pred. cov)');
set(gca,'ColorOrderIndex',1);
plot(t_vec,rmse_est,'--','DisplayName','RMSE (est. act.)');
plot(t_vec(2:end), rmse_pred, '--','DisplayName','RMSE (pred. act.)');
grid on;
xlabel('Time [sec]');
ylabel('Error [m]');
set(gca,'yscale','log');
utils.setPlotStyle(gca,{'widescreen'});

%% Repeat for Statistical Certainty
num_mc = 1000;
sse_pred = zeros(num_mc, num_time-1);
sse_est = zeros(num_mc, num_time);
sse_cov_pred = zeros(num_mc, num_time);
sse_cov_est = zeros(num_mc, num_time);
fprintf('Repeating tracker test for %d Monte Carlo trials...\n',num_mc);
num_trials_per_tic = 10;
num_tics_per_row=50;
for idx_mc = 1:num_mc
    if mod(idx_mc,num_trials_per_tic)==0
        fprintf('.');
    end
    if mod(idx_mc,num_trials_per_tic*num_tics_per_row)==0
        fprintf(' (%d/%d)\n',idx_mc,num_mc);
    end

    % Generate Measurements
    noise = L * randn(num_msmt, num_time);
    zeta = z + noise;
    
    % Initialize Track State
    x_pred = zeros(num_states,1);
    P_pred = zeros(num_states, num_states);
    
    % Initialize position with TDOA estimate from first measurement
    x_init = [0;50e3;5e3];
    epsilon = 100;
    x_pred(pos_idx) = tdoa.gdSoln(x_tdoa, zeta(:,1), C_roa, x_init,[],[],epsilon,[],[],[],ref_idx);
    P_pred(pos_idx,pos_idx) = tdoa.computeCRLB(x_tdoa,x_pred(pos_idx),C_roa,ref_idx,false,true);

    
    % Bound initial velocity uncertainty by assumed max velocity of 340 m/s
    % (Mach 1 at sea level)
    P_pred(vel_idx, vel_idx) = max_vel^2*eye(num_dims);
    
    % Step Through Time
    x_ekf_est = zeros(num_dims,num_time);
    x_ekf_pred = zeros(num_dims, num_time);
    for idx=1:num_time
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

        sse_cov_est(idx_mc, idx) = trace(P_est(pos_idx, pos_idx));
        sse_cov_pred(idx_mc, idx) = trace(P_pred(pos_idx, pos_idx));
    end

    err_pred = x_ekf_pred(:,1:end-1) - x_tgt_full(:,2:end);
    err_est = x_ekf_est - x_tgt_full;

    sse_pred(idx_mc,:) = sum(abs(err_pred).^2,1);
    sse_est(idx_mc,:) = sum(abs(err_est).^2,1);

end
fprintf('done.\n');

rmse_pred = sqrt(mean(sse_pred,1));
rmse_est = sqrt(mean(sse_est,1));

rmse_cov_est = sqrt(mean(sse_cov_est,1));
rmse_cov_pred = sqrt(mean(sse_cov_pred,1));

fig3=figure;
plot(t_vec,rmse_cov_est,'DisplayName','RMSE (est. cov.)');  
hold on;
plot(t_vec(2:end), rmse_cov_pred(1:end-1),'DisplayName','RMSE (pred. cov)');
set(gca,'ColorOrderIndex',1);
plot(t_vec,rmse_est,'--','DisplayName','RMSE (est. act.)');
plot(t_vec(2:end),rmse_pred,'--','DisplayName','RMSE (pred. act.)');
grid on;
xlabel('Time [sec]');
ylabel('Error [m]');
set(gca,'yscale','log');

legend('Location','NorthEast');
utils.setPlotStyle(gca,{'widescreen'});

%% Return Figure Handles
figs = [fig1, fig2, fig3];