function figs = book2_ex8_4()
% figs = book2_ex8_4()
%
% Executes Example 8.4 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% Constant-Turn Aircraft Tracking.
%
% An aircraft executes a sustained coordinated level turn (one full circle)
% tracked by four ground-based TDOA sensors.  Two EKF variants are
% compared:
%
%   CV model  - centripetal acceleration is unmodelled; large process noise
%               required, producing noisy estimates with visible lag.
%   CT model  - turn rate omega is a tracked state.  Small kinematic
%               process noise suffices; omega converges from 0 to the true
%               value within a few updates.
%
% INPUTS
%   none
%
% OUTPUTS
%   figs    [fig1, fig1b, fig2, fig3]  figure handles
%           fig1  – x-y (top-down) trajectory
%           fig1b – 3-D isometric trajectory (reveals altitude errors)
%           fig2  – 3-D position RMSE vs time
%           fig3  – estimated turn rate vs time
%
% Nicholas O'Donoughue
% June 2025

fprintf('Example 8.4...\n');

%% Target trajectory: constant-rate level turn ------------------------------
omega_true = pi / 60;               % rad/s  -> full circle in 120 s
v0         = 200;                   % m/s
alt        = 3e3;                   % m  (constant altitude)
R_turn     = v0 / omega_true;       % turn radius ~3820 m

t_inc  = 5.0;
t_max  = 2*pi / omega_true;         % 120 s  (one full circle)
t_vec  = (0 : t_inc : t_max - t_inc)';
num_time = numel(t_vec);

x_tgt_full = [ -R_turn*(1 - cos(omega_true*t_vec')); ...
                R_turn*      sin(omega_true*t_vec'); ...
                alt * ones(1, num_time) ];

%% TDOA sensors: rectangle surrounding the turn circle ----------------------
x_tdoa = [ 2e3, -9e3, -9e3,  2e3;
            5e3,  5e3, -5e3, -5e3;
             0.,   0.,   0.,   0.];
n_tdoa  = size(x_tdoa, 2);

sigma_toa = 100e-9;               % 10 ns  -> ~3 m RDOA noise
C_toa     = sigma_toa^2 * eye(n_tdoa);
C_roa     = utils.constants.c^2 * C_toa;
R         = utils.resampleCovMtx(C_roa, []);    % default ref sensor
L         = chol(R, 'lower');
num_msmt  = size(R, 1);
ref_idx   = [];

z    = tdoa.measurement(x_tdoa, x_tgt_full, ref_idx);
zeta = z + L * randn(num_msmt, num_time);

%% CV motion model ----------------------------------------------------------
% sigma_a = 1 m/s^2 (tuned for a non-maneuvering target).  The centripetal
% acceleration is ~10.5 m/s^2, so each 5-second step introduces a large
% model-mismatch error.
sigma_a_cv = 1;
motion_cv = tracker.makeMotionModel('cv', 3, sigma_a_cv^2);
F_cv = motion_cv.f_fun(t_inc);
Q_cv = motion_cv.q_fun(t_inc);
ss_cv = motion_cv.state_space;

[z_fun_cv, h_fun_cv] = tracker.makeMeasurementModel([], x_tdoa, [], [], ref_idx, [], ss_cv);

pos_idx_cv = ss_cv.pos_idx;
vel_idx_cv = ss_cv.vel_idx;
num_states_cv = ss_cv.num_states;

%% CT motion model ----------------------------------------------------------
% Same kinematic sigma_a; omega is treated as nearly-constant (small
% process_covar_omega keeps the estimate locked once it converges).
sigma_a_ct          = 1;
process_covar_omega = 1e-3;    % rad^2/s^3 spectral density

motion_ct = tracker.makeMotionModel('ct', 3, sigma_a_ct^2, process_covar_omega);
Q_ct = motion_ct.q_fun(t_inc);
f_ct = @(x) motion_ct.f_fun_ekf(x, t_inc);
g_ct = @(x) motion_ct.jacobian_fun(x, t_inc);

ss_ct = motion_ct.state_space;
[z_fun_ct, h_fun_ct] = tracker.makeMeasurementModel([], x_tdoa, [], [], ref_idx, [], ss_ct);

pos_idx_ct = ss_ct.pos_idx;
vel_idx_ct = ss_ct.vel_idx;
omega_idx  = ss_ct.omega_idx;
num_states_ct = ss_ct.num_states;

%% Shared initial offset (same for both trackers) ---------------------------
x_offset = [400; -400; 100];     % m  (position error at t=0)

% CV initial state
x0_cv = zeros(num_states_cv, 1);
P0_cv = zeros(num_states_cv, num_states_cv);
x0_cv(pos_idx_cv) = x_tgt_full(:, 1) + x_offset;
x0_cv(vel_idx_cv) = [0; v0; 0];
P0_cv(pos_idx_cv, pos_idx_cv) = (500^2) * eye(3);
P0_cv(vel_idx_cv, vel_idx_cv) = (100^2) * eye(3);

% CT initial state  (omega unknown -> start at 0 with generous uncertainty)
sigma_omega_init = 0.1;          % rad/s
x0_ct = zeros(num_states_ct, 1);
P0_ct = zeros(num_states_ct, num_states_ct);
x0_ct(pos_idx_ct) = x_tgt_full(:, 1) + x_offset;
x0_ct(vel_idx_ct) = [0; v0; 0];
x0_ct(omega_idx)  = 0;
P0_ct(pos_idx_ct, pos_idx_ct) = (500^2) * eye(3);
P0_ct(vel_idx_ct, vel_idx_ct) = (100^2) * eye(3);
P0_ct(omega_idx, omega_idx)   = sigma_omega_init^2;

%% EKF loop ------------------------------------------------------------------
x_cv_est  = zeros(3, num_time);
x_ct_est  = zeros(3, num_time);
omega_est = zeros(1, num_time);

x_cv_pred  = x0_cv;  P_cv_pred  = P0_cv;
x_ct_pred  = x0_ct;  P_ct_pred  = P0_ct;

fprintf('Iterating through EKF tracker time steps...');
for idx = 1 : num_time
    fprintf('.');
    this_zeta = zeta(:, idx);

    % EKF Update
    [x_cv_k, P_cv_k] = tracker.ekfUpdate(x_cv_pred, P_cv_pred, this_zeta, R, z_fun_cv, h_fun_cv);
    [x_ct_k, P_ct_k] = tracker.ekfUpdate(x_ct_pred, P_ct_pred, this_zeta, R, z_fun_ct, h_fun_ct);

    % Store estimated positions and turn rate
    x_cv_est(:, idx)  = x_cv_k(pos_idx_cv);
    x_ct_est(:, idx)  = x_ct_k(pos_idx_ct);
    omega_est(idx)    = x_ct_k(omega_idx);

    % CV: linear KF predict
    [x_cv_pred, P_cv_pred] = tracker.kfPredict(x_cv_k, P_cv_k, Q_cv, F_cv);

    % CT: nonlinear EKF predict using exact CT transition + analytical Jacobian
    [x_ct_pred, P_ct_pred] = tracker.ekfPredict(x_ct_k, P_ct_k, Q_ct, ...
                                                f_ct, g_ct);
end
fprintf('done.\n');

%% Figure 1: x-y (top-down) trajectory -------------------------------------
fig1 = figure;
plot(x_tgt_full(1,:)/1e3, x_tgt_full(2,:)/1e3, 'k-', 'LineWidth', 2, ...
     'DisplayName', 'True trajectory');
hold on;
plot(x_cv_est(1,:)/1e3, x_cv_est(2,:)/1e3, '--', 'DisplayName', 'CV model');
plot(x_ct_est(1,:)/1e3, x_ct_est(2,:)/1e3, '-.', 'DisplayName', 'CT model');
plot(x_tdoa(1,:)/1e3, x_tdoa(2,:)/1e3, '^', 'MarkerFaceColor', 'auto', ...
     'DisplayName', 'Sensors');
xlabel('x [km]');
ylabel('y [km]');
axis equal;
legend('Location', 'NorthEast');
grid on;
% title('Example 8.4: Constant-Turn Aircraft – x-y Trajectory');
utils.setPlotStyle(gca, {'widescreen'});

%% Figure 1b: 3-D isometric trajectory (reveals altitude errors) ------------
fig1b = figure;
ax = axes('Parent', fig1b);
plot3(ax, x_tgt_full(1,:)/1e3, x_tgt_full(2,:)/1e3, x_tgt_full(3,:)/1e3, ...
      'k-', 'LineWidth', 2, 'DisplayName', 'True trajectory');
hold(ax, 'on');
plot3(ax, x_cv_est(1,:)/1e3, x_cv_est(2,:)/1e3, x_cv_est(3,:)/1e3, ...
      '--', 'DisplayName', 'CV model');
plot3(ax, x_ct_est(1,:)/1e3, x_ct_est(2,:)/1e3, x_ct_est(3,:)/1e3, ...
      '-.', 'DisplayName', 'CT model');
scatter3(ax, x_tdoa(1,:)/1e3, x_tdoa(2,:)/1e3, x_tdoa(3,:)/1e3, ...
         60, 'filled', '^', 'DisplayName', 'Sensors');
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
view(ax, 30, 45);   % isometric-like viewing angle (azimuth=30, elevation=45)
legend(ax, 'Location', 'NorthEast');
grid(ax, 'on');
% title('Example 8.4: Constant-Turn Aircraft – 3-D Trajectory');
utils.setPlotStyle(ax, {'widescreen'});

%% Figure 2: 3-D position RMSE vs time -------------------------------------
rmse_cv = sqrt(sum((x_cv_est - x_tgt_full).^2, 1));
rmse_ct = sqrt(sum((x_ct_est - x_tgt_full).^2, 1));

fig2 = figure;
semilogy(t_vec, rmse_cv/1e3, 'DisplayName', 'CV model');
hold on;
semilogy(t_vec, rmse_ct/1e3, 'DisplayName', 'CT model');
xlabel('Time [s]');
ylabel('3-D RMSE [km]');
legend('Location', 'NorthEast');
grid on;
% title('Example 8.4: Position RMSE vs Time');
utils.setPlotStyle(gca, {'widescreen'});

%% Figure 3: estimated turn rate vs time -----------------------------------
fig3 = figure;
yline(omega_true * 180/pi, 'k-', 'LineWidth', 2, 'DisplayName', 'True \omega');
hold on;
plot(t_vec, omega_est * 180/pi, 'DisplayName', 'CT estimate');
xlabel('Time [s]');
ylabel('Turn rate [deg/s]');
legend('Location', 'SouthEast');
grid on;
% title('Example 8.4: CT Model – Estimated Turn Rate');
utils.setPlotStyle(gca, {'widescreen'});

figs = [fig1, fig1b, fig2, fig3];
