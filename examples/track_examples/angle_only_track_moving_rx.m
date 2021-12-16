% Example of an angle-only track from a moving receiver
addpath('C:/Users/nodonoug/Documents/GitLab/emitter-detection-book/');

%% Target position
x_init = [100,90,10]'*1e3;
vt = [10,0,0]'; 

%% Sensor position
xs_init = [0,0,0; 0, 10,20; 30, 30,20]*1e3; % m
vs = [0,40,0]'; % m/s
n_sensor=size(xs_init,2);

%% Sensor performance -- assume 2D DOA
do_good_init = true;
good_init_rmse = 10e3;

% test_case = 'df_state_space';
% test_case = 'df_msmt_space';
% test_case = 'df_singlerx_msmt_space';
test_case = 'tdoa_msmt_space';
switch lower(test_case)        
    case 'df_state_space'
%         2D DF -- state space
        az_err_deg = 1;
        el_err_deg = 1;
        df_cov = diag(kron([az_err_deg, el_err_deg]*pi/180,ones(1,n_sensor))).^2;
        noise_premult = chol(df_cov,'upper');
        n_msmt = size(df_cov,1);
        msmt_model = @(xs, xt) triang.lsSoln(xs,triang.measurement(xs,xt)+noise_premult*randn(n_msmt,1),df_cov,mean([xs,xt],2));
        msmt_jacob_model = @(xs,xt) eye(size(xs,1)); 

        msmt_cov_model = @(xs,xt) triang.computeCRLB(xs, xt, df_cov);

        solver = @(xs,z,x_init) z;
    case 'df_msmt_space'
        msmt_model = @(xs,xt) triang.measurement(xs,xt,true);
        az_err_deg = 1;
        el_err_deg = 1;
        df_cov = diag(kron([az_err_deg, el_err_deg]*pi/180,ones(1,n_sensor))).^2;
        noise_premult = chol(df_cov,'upper');
        n_msmt = size(df_cov,1);
        solver = @(xs,z,x_init) triang.lsSoln(xs,z,df_cov,x_init);

        msmt_jacob_model = @(xs,xt) triang.jacobian(xs,xt,true);     
        msmt_cov_model = @(xs,xt) triang.computeCRLB(xs, xt, df_cov);
    case 'df_singlerx_msmt_space'
        xs_init = xs_init(:,1);
        xs_init(2) = -400e3;
        xs_init(3) = 500e3;
        vs = [0,200,0]';
        n_sensor = 1;
        
        msmt_model = @(xs,xt) triang.measurement(xs,xt,true);
        az_err_deg = 1;
        el_err_deg = 1;
        df_cov = diag(kron([az_err_deg, el_err_deg]*pi/180,ones(1,n_sensor))).^2;
        noise_premult = chol(df_cov,'upper');
        n_msmt = size(df_cov,1);
        solver = @(xs,z,x_init) triang.lsSoln(xs,z,df_cov,x_init);

        msmt_jacob_model = @(xs,xt) triang.jacobian(xs,xt,true);     
        msmt_cov_model = @(xs,xt) triang.computeCRLB(xs, xt, df_cov);
    case 'tdoa_msmt_space'
        % TDOA
        time_err = 1e-7;
        roa_cov = (time_err*utils.constants.c)^2*eye(n_sensor);
        tdoa_ref_idx = 1;
        rdoa_cov = utils.resampleCovMtx(roa_cov,tdoa_ref_idx);
        noise_premult = chol(rdoa_cov,'upper');
        n_msmt = size(rdoa_cov,2);
        % 
        msmt_model = @(xs,xt) tdoa.measurement(xs,xt,tdoa_ref_idx);
        msmt_jacob_model = @(xs,xt) tdoa.jacobian(xs,xt,tdoa_ref_idx);
        msmt_cov_model = @(xs,xt) tdoa.computeCRLB(xs, xt, roa_cov,tdoa_ref_idx,false,true);
        solver = @(xs,z,x_init) tdoa.lsSoln(xs,z,roa_cov,x_init,[],[],[],[],tdoa_ref_idx);
end


%% Tracker Assumptions
% We're building this tracker in measurement space, rather than state
% space, so use the direct measurement error covariance, rather than the
% position estimation CRLB from that measurement.  This is necessary if
% individual measurements lack the power for generating a geolocation
% result.

update_rate = 60; % seconds
heading = 0;
latency = 0; % processing latency
maxG=1; % assumed maneuverability
num_states = 2; % pos/vel/accel
num_dims = 3;

% Make State Matrices
F = tracker.makeTransitionMatrix(update_rate, num_states, num_dims);
    % States are: x, vx, ax, y, vy, zy, z, vx, az

% Directly measure x, y, and z
H = kron(eye(num_dims),[1,zeros(1,num_states-1)]);

target_type = 'lateral-only';
Q = tracker.makeCAProcessNoise(maxG, num_states, heading, 0, update_rate, target_type);
%R = df_cov;
R = rdoa_cov;

%% Scenario Assumptions
maxT = 60*60; % 60 minutes, represented in seconds
t = 0:update_rate:maxT;
num_frames = numel(t);

x_true_full = zeros(num_dims,num_frames);
x_msmt_full = zeros(num_dims,num_frames);
x_est_full = zeros(num_dims,num_frames);
x_pred_full = zeros(num_dims,num_frames);
x_sensor_full = zeros(num_dims,n_sensor,num_frames);
P_est_full = zeros(num_dims,num_dims,num_frames);
P_pred_full = zeros(num_dims,num_dims,num_frames);

% Initialize Scene and Tracker State
xt = x_init;
xs = xs_init;
pos_mask = logical(kron(ones(1,num_dims),[1,zeros(1,num_states-1)])); % binary mapping to the x, y, and z entries in the state matrix

% Good initialization (e.g. from off-board early warning)
if do_good_init
    % Define covariance
    Ppos = good_init_rmse^2/num_dims * eye(num_dims);
    P_pred = kron(Ppos,[1,0;0,10*eye(num_states-1)]);
    
    % Random initialization
    x_pred = zeros(num_states*num_dims,1);
    x_pred(pos_mask) = xt + sqrt(Ppos)*randn(num_dims,1);
else
    x_pred = kron([10, 10, 10],[1,zeros(1,num_states-1)])'*1e3;
    P_pred = diag(kron([100e3,100e3,5e3],[1,.5*ones(1,num_states-1)])).^2;
end

% Assume level flight --- manually set P(6,6) to zero.
P_pred(:,end) = 0;
P_pred(end,:) = 0;
P_pred(end,end) = .01;

%% Handle Warnings
warning('off','MATLAB:nearlySingularMatrix')

%% Run Scenario
for ii = 1:num_frames
    % Generate Noisy Measurement
    switch lower(test_case)
        case {'df_msmt_space','df_singlerx_msmt_space'}
            R = df_cov;
            
            z_opt = msmt_model(xs,xt);
            n = noise_premult * randn(n_msmt,1);
            z = z_opt + n;
            x_msmt_full(:,ii) = solver(xs,z,x_pred(pos_mask));
            
        case 'df_state_space'
            R = msmt_cov_model(xs,xt);

            % Error already in measurement model; no need to add it
            z = msmt_model(xs,xt);
            
            % Measurement in state space; no need to "solve" for x
            x_msmt_full(:,ii) = z;
            
        case 'tdoa_msmt_space'
            R = rdoa_cov;
           
            z_opt = msmt_model(xs,xt);
            n = noise_premult * randn(n_msmt,1);
            z = z_opt + n;
            x_msmt_full(:,ii) = solver(xs,z,x_pred(pos_mask));
    end

    
   % Initialize tracker from DF geolocation, if first frame
   if ii==1
       if ~do_good_init && ~any(isnan(x_msmt_full(:,ii))) % Insufficient info to re-initialize
           x_pred(pos_mask) = x_msmt_full(:,ii);
           P_pred(pos_mask,pos_mask) = msmt_cov_model(xs,x_pred(pos_mask));
       end
       fprintf('Initializing track with single update geolocation...\n');
       fprintf('Expected error: %.2f km\n',sqrt(trace(P_pred(pos_mask,pos_mask)))/1e3);
       fprintf('Actual error: %.2f km\n',sqrt(sum(abs(x_pred(pos_mask)-xt).^2))/1e3);
       fprintf('Update, Predicted Err, Measured Err, Post-Update Err\n');
   end
   
   % Store prediction from prior update and current truth
   x_pred_full(:,ii) = x_pred(pos_mask);
   P_pred_full(:,:,ii) = P_pred(pos_mask,pos_mask);
   x_true_full(:,ii) = xt;
   x_sensor_full(:,:,ii) = xs;

   if ii==1, continue, end
   
   % Update EKF Measurement Matrix
    this_H = @(x) msmt_model(xs,x(pos_mask));
    this_H_jacob = @(x) kron(msmt_jacob_model(xs,x(pos_mask)), [1,zeros(1,num_states-1)]')';
   
   % Run KF Step
   [x_pred, x_est, P_pred, P_est] = tracker.runKalmanStep(x_pred, P_pred, z, F, this_H, Q, R, 'kf', 'ekf', [], this_H_jacob);
%    [x_pred, x_est, P_pred, P_est] = tracker.runKalmanStep(x_pred, P_pred, z, F, H, Q, R, 'kf', 'kf');
   
   % Store estimated position
   x_est_full(:,ii) = x_est(pos_mask);
   P_est_full(:,:,ii) = P_est(pos_mask,pos_mask);
   
   % Report on errors
   err_msmt = x_msmt_full(:,ii) - xt;
   err_pred = x_pred_full(:,ii) - xt;
   err_est = x_est(pos_mask) - xt;
   
   fprintf('%6d,       %.2f km,      %.2f km,         %.2f km\n', ii, ...
            sqrt(sum(abs(err_pred).^2))/1e3, ...
            sqrt(sum(abs(err_msmt).^2))/1e3, ...
            sqrt(sum(abs(err_est).^2))/1e3);
   % Update positions
   xt = xt + update_rate*vt;
   xs = xs + update_rate*vs;
end

warning('on','MATLAB:nearlySingularMatrix')

%% Post Process
rmse_pred = zeros(1,num_frames);
rmse_est = zeros(1,num_frames);
rmse_P_pred = zeros(1,num_frames);
rmse_P_est = zeros(1,num_frames);
rmse_msmt = zeros(1,num_frames);

rmse_est_xy = zeros(1,num_frames);

for ii=1:num_frames
    err_pred = x_pred_full(:,ii) - x_true_full(:,ii);
    err_est = x_est_full(:,ii) - x_true_full(:,ii);
    
    rmse_pred(ii) = norm(err_pred)/1e3;
    rmse_est(ii) = norm(err_est)/1e3;
    rmse_est_xy(ii) = norm(err_est(1:2))/1e3;
    
    rmse_P_pred(ii) = sqrt(sum(diag(P_pred_full(:,:,ii))))/1e3;
    rmse_P_est(ii) = sqrt(sum(diag(P_est_full(:,:,ii))))/1e3;
    
    err_msmt = x_msmt_full(:,ii) - x_true_full(:,ii);
    rmse_msmt(ii) = norm(err_msmt)/1e3;
end

%% Plots
figure;
hdl=plot(t,rmse_P_pred);
hold on;
plot(t,rmse_pred,'--','Color',hdl.Color);
hdl=plot(t,rmse_P_est);
plot(t,rmse_est,'--','Color',hdl.Color);
plot(t,rmse_est_xy,'-.');
plot(t,rmse_msmt)
legend('Prediction Error (theo.)','Prediction Error (actual)',...
       'Estimation Error (theo.)','Estimation Error (actual)',...
       'Estimation Error (x-y only)','Msmt Error');
xlabel('Time [s]');
ylabel('Error [km]');
set(gca,'yscale','log')
ylim([.1, 1e3]);
grid on;


figure;
plot3(squeeze(x_true_full(1,:,:))/1e3,squeeze(x_true_full(2,:,:))/1e3,...
      squeeze(x_true_full(3,:,:))/1e3,'DisplayName','True');
hold on;
plot3(x_sensor_full(1,:)/1e3,x_sensor_full(2,:)/1e3,x_sensor_full(3,:)/1e3,'DisplayName','Sensor');
plot3(x_pred_full(1,:)/1e3,x_pred_full(2,:)/1e3,x_pred_full(3,:)/1e3,'DisplayName','Prediction');
plot3(x_est_full(1,:)/1e3, x_est_full(2,:)/1e3, x_est_full(3,:)/1e3, 'DisplayName','Estimate');
plot3(x_msmt_full(1,:)/1e3,x_msmt_full(2,:)/1e3,x_msmt_full(3,:)/1e3,'DisplayName','Msmt Solution (seeded by tracker)');
grid on;
xlabel('Cross-Track [km]');
ylabel('Along-Track [km]');
zlabel('Altitude [km]');
xlim([0,2]*max(x_init(1))/1e3);
ylim([-1,1]*max(max(abs(x_sensor_full(2,:,:))))/1e3);
zlim([0, max(xs_init(3,:))]/1e3);

figure;
plot3(squeeze(x_true_full(1,:,:))/1e3,squeeze(x_true_full(2,:,:))/1e3,...
      squeeze(x_true_full(3,:,:))/1e3,'DisplayName','True');
hold on;
plot3(x_sensor_full(1,:)/1e3,x_sensor_full(2,:)/1e3,x_sensor_full(3,:)/1e3,'DisplayName','Sensor');
plot3(x_pred_full(1,:)/1e3,x_pred_full(2,:)/1e3,x_pred_full(3,:)/1e3,'DisplayName','Prediction');
plot3(x_est_full(1,:)/1e3, x_est_full(2,:)/1e3, x_est_full(3,:)/1e3, 'DisplayName','Estimate');
plot3(x_msmt_full(1,:)/1e3,x_msmt_full(2,:)/1e3,x_msmt_full(3,:)/1e3,'DisplayName','Msmt Solution (seeded by tracker)');
grid on;
xlabel('Cross-Track [km]');
ylabel('Along-Track [km]');
zlabel('Altitude [km]');
xlim([0,2]*max(x_init(1))/1e3);
ylim([-1,1]*max(max(abs(x_sensor_full(2,:,:))))/1e3);
zlim([0, 2*max(max(x_true_full(3,:,:)))]/1e3);