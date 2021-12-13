% Example of an angle-only track from a moving receiver


%% Target position
x_init = [100,90,10]'*1e3;
vt = [10,0,0]'; 

%% Sensor position
xs_init = [0,0,0,0; 0, 10,20,30; 30, 30,20,20]*1e3; % m
vs = [0,20,0]'; % m/s
n_sensor=size(xs_init,2);

%% Sensor performance -- assume 2D DOA

% 2D DF -- state space
% az_err_deg = 1;
% el_err_deg = 1;
% df_cov = diag(kron([az_err_deg, el_err_deg]*pi/180,ones(1,n_sensor))).^2;
% noise_premult = chol(df_cov,'upper');
% n_msmt = size(df_cov,1);
% msmt_model = @(xs, xt) triang.lsSoln(xs,triang.measurement(xs,xt)+noise_premult*randn(n_msmt,1),df_cov,mean([xs,xt],2));
% msmt_jacob_model = @(xs,xt) eye(size(xs,1)); 
% 
% msmt_cov_model = @(xs,xt) triang.computeCRLB(xs, xt, df_cov);
% 
% solver = @(xs,z,x_init) z;

% % TDOA
% time_err = 1e-9;
% tdoa_cov = (time_err*utils.constants.c)^2*(1+eye(n_sensor-1));
% noise_premult = chol(tdoa_cov,'upper');
% tdoa_ref_idx = 1;
% n_msmt = size(tdoa_cov,2);
% % 
% msmt_model = @(xs,xt) tdoa.measurement(xs,xt,tdoa_ref_idx);
% msmt_jacob_model = @(xs,xt) tdoa.jacobian(xs,xt,tdoa_ref_idx);
% solver = @(xs,z,x_init) tdoa.lsSoln(xs,z,tdoa_cov,x_init,[],[],[],[],tdoa_ref_idx);

msmt_model = @(xs,xt) triang.measurement(xs,xt,true);
az_err_deg = .01;
el_err_deg = .01;
df_cov = diag(kron([az_err_deg, el_err_deg]*pi/180,ones(1,n_sensor))).^2;
noise_premult = chol(df_cov,'upper');

solver = @(xs,z,x_init) triang.lsSoln(xs,z,df_cov,x_init);

% msmt_err_model = @(xs,xt,C) triang.computeCRLB(xs,xt,C);
msmt_jacob_model = @(xs,xt) triang.jacobian(xs,xt,true);

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

H = cat(2,eye(num_dims), zeros(num_dims,(num_states-1)*num_dims));

target_type = 'lateral-only';
Q = tracker.makeCAProcessNoise(maxG, num_states, heading, 0, update_rate, target_type);
R = df_cov;
%R = tdoa_cov;

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
x_pred = [10, 10, 10]'*1e3;
P_pred = diag([100e3,100e3,5e3]).^2;
if num_states > 1
    x_pred = cat(1,x_pred,zeros(num_dims*(num_states-1),1));
    P_pred = [P_pred, zeros(num_dims,num_dims*(num_states-1));
              zeros(num_dims*(num_states-1),num_dims), 200*eye(num_dims*(num_states-1))];
end

%% Run Scenario
for ii = 1:num_frames
   % Generate Noisy Measurement
   z_opt = msmt_model(xs,xt);
%    R = msmt_cov_model(xs,xt);
   n = noise_premult * randn(n_msmt,1);
   z = z_opt + n;
   x_msmt_full(:,ii) = solver(xs,z,x_pred(1:num_dims));
    
   % Initialize tracker from DF geolocation, if first frame
   if ii==1
       x_pred(1:num_dims) = x_msmt_full(:,ii);
       P_pred(1:num_dims,1:num_dims) = triang.computeCRLB(xs,x_pred(1:num_dims),df_cov);
       fprintf('Initializing track with single update geolocation...\n');
       fprintf('Expected error: %.2f km\n',sqrt(trace(P_pred(1:num_dims,1:num_dims)))/1e3);
       fprintf('Actual error: %.2f km\n',sqrt(sum(abs(x_pred(1:num_dims)-xt).^2))/1e3);
   end
   
   % Store prediction from prior update and current truth
   x_pred_full(:,ii) = x_pred(1:num_dims);
   P_pred_full(:,:,ii) = P_pred(1:num_dims,1:num_dims);
   x_true_full(:,ii) = xt;
   x_sensor_full(:,:,ii) = xs;

   if ii==1, continue, end
   
   % Update EKF Measurement Matrix
%    this_H = @(x) msmt_model(xs,x(1:num_dims));
%    this_H_jacob = @(x) cat(1,msmt_jacob_model(xs,x(1:num_dims)), zeros(num_dims*(num_states-1),n_msmt))';
   
   % Run KF Step
   [x_pred, x_est, P_pred, P_est] = tracker.runKalmanStep(x_pred, P_pred, z, F, this_H, Q, R, 'kf', 'ekf', [], this_H_jacob);
%    [x_pred, x_est, P_pred, P_est] = tracker.runKalmanStep(x_pred, P_pred, z, F, H, Q, R, 'kf', 'kf');
   
   % Store estimated position
   x_est_full(:,ii) = x_est(1:num_dims);
   P_est_full(:,:,ii) = P_est(1:num_dims,1:num_dims);
   
   % Update positions
   xt = xt + update_rate*vt;
   xs = xs + update_rate*vs;
end

%% Post Process
rmse_pred = zeros(1,num_frames);
rmse_est = zeros(1,num_frames);
rmse_P_pred = zeros(1,num_frames);
rmse_P_est = zeros(1,num_frames);
rmse_msmt = zeros(1,num_frames);

for ii=1:num_frames
    err_pred = x_pred_full(:,ii) - x_true_full(:,ii);
    err_est = x_est_full(:,ii) - x_true_full(:,ii);
    
    rmse_pred(ii) = norm(err_pred)/1e3;
    rmse_est(ii) = norm(err_est)/1e3;
    
    rmse_P_pred(ii) = sqrt(sum(diag(P_pred)))/1e3;
    rmse_P_est(ii) = sqrt(sum(diag(P_est)))/1e3;
    
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
plot(t,rmse_msmt)
legend('Prediction Error (theo.)','Prediction Error (actual)',...
       'Estimation Error (theo.)','Estimation Error (actual)','Msmt Error');
xlabel('Time [s]');
ylabel('Error [km]');
set(gca,'yscale','log')

figure;
plot3(squeeze(x_true_full(1,:,:)),squeeze(x_true_full(2,:,:)),...
      squeeze(x_true_full(3,:,:)),'DisplayName','True');
hold on;
plot3(x_sensor_full(1,:),x_sensor_full(2,:),x_sensor_full(3,:),'DisplayName','Sensor');
plot3(x_pred_full(1,:),x_pred_full(2,:),x_pred_full(3,:),'DisplayName','Prediction');
plot3(x_est_full(1,:),x_est_full(2,:),x_est_full(3,:),'DisplayName','Estimate');
grid on;