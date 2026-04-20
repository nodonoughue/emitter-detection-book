function figs = book2_ex8_3()
% figs = book2_ex8_3()
%
% Executes Example 8.3 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% Ballistic Trajectory Tracking.
%
% A projectile is launched and tracked using four ground-based TDOA
% sensors.  Two EKF variants are compared:
%
%   CV model    - gravity is unmodelled; requires large process noise to
%                 follow the curved descent.
%   Ballistic   - gravity is an explicit deterministic forcing term on the
%                 state mean.  Small kinematic process noise suffices.
%
% INPUTS
%   none
%
% OUTPUTS
%   figs    [fig1, fig2]  figure handles
%
% Nicholas O'Donoughue
% June 2025

fprintf('Example 8.3...\n');

%% Target trajectory: ballistic, 3-D ----------------------------------------
g       = -9.80665;                  % m/s^2  (downward, z-axis)
t_inc   = 1.0;                       % s  between updates
v_init  = [100; 80; 200];            % initial velocity [m/s]  east, north, up
x_init  = [0; 0; 0];                 % launch point (ground level)

t_flight = -2 * v_init(3) / g;      % total flight time ~40.8 s
t_vec    = (0 : t_inc : t_flight - t_inc)';
num_time = numel(t_vec);

gravity_vec  = [0; 0; g];
x_tgt_full   = x_init + v_init .* t_vec' + 0.5 * gravity_vec .* (t_vec.^2)';
% size: 3 x num_time

%% TDOA sensors: cross pattern at 3 km --------------------------------------
x_tdoa = [ 3e3,  0., -3e3,  0.;
            0.,  3e3,  0., -3e3;
           10.,  10.,  10.,  10.];
n_tdoa  = size(x_tdoa, 2);

ref_idx   = 1;
sigma_toa = 300e-9;              % 300 ns  -> ~90 m RDOA noise
C_toa     = sigma_toa^2 * eye(n_tdoa);
C_roa     = utils.constants.c^2 * C_toa;
R         = utils.resampleCovMtx(C_roa, ref_idx);
L         = chol(R, 'lower');
num_msmt  = size(R, 1);

z    = tdoa.measurement(x_tdoa, x_tgt_full, ref_idx);
zeta = z + L * randn(num_msmt, num_time);

%% Motion models -------------------------------------------------------------
% CV: large process noise to compensate for unmodelled gravity (~50 m/s^2)
q_a_cv = 50;
motion_cv = tracker.makeMotionModel('cv', 3, q_a_cv^2);
F_cv = motion_cv.f_fun(t_inc);
Q_cv = motion_cv.q_fun(t_inc);

% Ballistic: same CV dynamics; gravity is applied as a mean offset after
% kfPredict.  Small kinematic noise suffices because gravity is explicit.
q_a_bal = 3;
motion_bal = tracker.makeMotionModel('ballistic', 3, q_a_bal^2);
F_bal = motion_bal.f_fun(t_inc);
Q_bal = motion_bal.q_fun(t_inc);
b_bal = motion_bal.b_fun(t_inc);

ss_cv = motion_cv.state_space;

% Both share the same TDOA measurement model (same state-space layout)
msmt = tracker.makeMeasurementModel([], x_tdoa, [], [], ref_idx, []);

pos_idx = ss_cv.pos_idx;
vel_idx = ss_cv.vel_idx;
num_states = ss_cv.num_states;

%% Shared initial conditions (intentionally perturbed) ----------------------
x0     = zeros(num_states, 1);
P0     = zeros(num_states, num_states);
x0(pos_idx) = [300; -200; 200];   % perturbed from the true launch point
x0(vel_idx) = [50;  50;   100];   % wrong initial velocity
P0(pos_idx, pos_idx) = (500^2) * eye(3);
P0(vel_idx, vel_idx) = (300^2) * eye(3);

%% EKF loop ------------------------------------------------------------------
x_cv_est  = zeros(3, num_time);
x_bal_est = zeros(3, num_time);

x_cv_pred  = x0;  P_cv_pred  = P0;
x_bal_pred = x0;  P_bal_pred = P0;

fprintf('Iterating through EKF tracker time steps...');
for idx = 1 : num_time
    fprintf('.');
    this_zeta = zeta(:, idx);

    % EKF Update
    [x_cv_est_k,  P_cv_est]  = tracker.ekfUpdate(x_cv_pred,  P_cv_pred,  this_zeta, R, msmt.z_fun_raw, msmt.h_fun_raw);
    [x_bal_est_k, P_bal_est] = tracker.ekfUpdate(x_bal_pred, P_bal_pred, this_zeta, R, msmt.z_fun_raw, msmt.h_fun_raw);

    % Store estimated positions
    x_cv_est(:, idx)  = x_cv_est_k(pos_idx);
    x_bal_est(:, idx) = x_bal_est_k(pos_idx);

    % KF Predict (linear step; ballistic model passes gravity as control input u)
    [x_cv_pred,  P_cv_pred]  = tracker.kfPredict(x_cv_est_k,  P_cv_est,  Q_cv,  F_cv);
    [x_bal_pred, P_bal_pred] = tracker.kfPredict(x_bal_est_k, P_bal_est, Q_bal, F_bal, b_bal);
end
fprintf('done.\n');

%% Figure 1: range-altitude profile -----------------------------------------
rng_tgt = hypot(x_tgt_full(1,:), x_tgt_full(2,:));
rng_cv  = hypot(x_cv_est(1,:),  x_cv_est(2,:));
rng_bal = hypot(x_bal_est(1,:), x_bal_est(2,:));

fig1 = figure;
plot(rng_tgt/1e3, x_tgt_full(3,:)/1e3, 'k-', 'LineWidth', 2, ...
     'DisplayName', 'True trajectory');
hold on;
plot(rng_cv/1e3,  x_cv_est(3,:)/1e3,  '--', 'DisplayName', 'CV model');
plot(rng_bal/1e3, x_bal_est(3,:)/1e3, '-.', 'DisplayName', 'Ballistic model');
xlabel('Horizontal range [km]');
ylabel('Altitude [km]');
legend('Location', 'NorthEast');
grid on;
% title('Example 8.3: Ballistic Trajectory – Range/Altitude Profile');
utils.setPlotStyle(gca, {'widescreen'});

%% Figure 2: 3-D position RMSE vs time --------------------------------------
rmse_cv  = sqrt(sum((x_cv_est  - x_tgt_full).^2, 1));
rmse_bal = sqrt(sum((x_bal_est - x_tgt_full).^2, 1));

fig2 = figure;
semilogy(t_vec, rmse_cv/1e3,  'DisplayName', 'CV model');
hold on;
semilogy(t_vec, rmse_bal/1e3, 'DisplayName', 'Ballistic model');
xlabel('Time [s]');
ylabel('3-D RMSE [km]');
legend('Location', 'NorthEast');
grid on;
% title('Example 8.3: Position RMSE vs Time');
utils.setPlotStyle(gca, {'widescreen'});

figs = [fig1, fig2];
