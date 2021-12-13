%% Setup Tracker Parameters

fprintf('Running steady-state tracker example...\n');

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

% Measurement Error
R = diag([10e3, 10e3, 20e3].^2); % ENU
H = H_pos;

% Predict track error using horizontal divert model (steady-state,
% level flight)
% Rng-only PCL, No AOA, No Doppler
[~, Pe] = tracker.steadystateError(F,H,Q,R);

%% Print Output
fprintf('Steady-State Position Prediction Error:\n');
pos_idx = 1:num_dims;
disp(Pe(pos_idx,pos_idx));
if num_states > 1
    fprintf('Steady-State Velocity Prediction Error:\n');
    vel_idx = num_dims + (1:num_dims);
    disp(Pe(vel_idx, vel_idx));
end
if num_states > 2
    fprintf('Steady-State Acceleration Prediction Error:\n');
    acc_idx = num_dims + (1:num_dims);
    disp(Pe(acc_idx, acc_idx));
end

