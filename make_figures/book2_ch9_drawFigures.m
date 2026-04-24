% Draw Figures - Chapter 9
%
% This script generates all of the figures that appear in
% Chapter 9 of the textbook.
%
% Nicholas O'Donoughue
% June 2025

% Clear Figures
close all;

% Set up directory and filename for figures
dirNm = fullfile(pwd,'figures','practical_geo');
if ~exist(dirNm,'dir')
    mkdir(dirNm);
end
prefix = fullfile(dirNm,'fig9_');

% Initialize Plot Preference
doFullColor=true;
utils.initPlotSettings;
colors=get(groot,'DefaultAxesColorOrder');

% Reset the random number generator, to ensure reproducibility
rng('default');

addpath('examples');

%% Figure 9.2 - Track, Prediction, and Acceptance Gate
%
% One 4-state CV track in 2D TDOA space; predict one step forward and
% show the acceptance gate in both x/y (position) and zeta (TDOA) space.
% Ten random measurement states are scattered around the prediction and
% classified as valid (inside gate) or invalid (outside gate).

fprintf('Generating Figure 9.2...\n');

% Motion model: CV, 2D, process noise 100 m^2/s^4 per dim
mm_f2  = tracker.makeMotionModel('cv', 2, 100*eye(2));
ss_f2  = mm_f2.state_space;
dt_f2  = 0.5;    % [s] time step

% Track states at t = 0, 0.5, 1.0, 1.5  (columns: [px; py; vx; vy])
sv_f2 = 1e3 * [0.00, 0.50, 1.00, 1.25;   % px
                0.00, 0.50, 0.50, 1.00;   % py
                1.00, 1.00, 0.50, 0.50;   % vx
                1.00, 0.00, 1.00, 0.50];  % vy
P_last_f2 = diag([2.0, 0.75, 0.4, 0.4]) * 1e4;

trk_f2 = tracker.makeTrack(tracker.makeState(ss_f2, 0*dt_f2, sv_f2(:,1)), '0');
for kk = 2:3
    trk_f2 = tracker.appendTrack(trk_f2, tracker.makeState(ss_f2, (kk-1)*dt_f2, sv_f2(:,kk)));
end
trk_f2 = tracker.appendTrack(trk_f2, tracker.makeState(ss_f2, 3*dt_f2, sv_f2(:,4), P_last_f2));

% Predict to t = 2.0 s
s_pred_f2 = tracker.predictState(tracker.currState(trk_f2), 4*dt_f2, mm_f2);

% TDOA sensor array (3 sensors, 2D)
x_tdoa_f2 = [-1e3, -1e3, -1e3;   % x [m]
               0,  500, 1e3];     % y [m]
ref_f2    = 1;
C_roa_f2  = 100 * eye(3);         % range variance [m^2]
R_f2      = utils.resampleCovMtx(C_roa_f2, ref_f2);    % 2x2

msmt_f2 = tracker.makeMeasurementModel([], x_tdoa_f2, [], [], ref_f2, [], R_f2);

% Predicted measurement and innovation covariance
z_pred_f2 = msmt_f2.z_fun(s_pred_f2);
H_f2      = msmt_f2.h_fun(s_pred_f2);
if ndims(H_f2) == 3; H_f2 = H_f2(:,:,1)'; end
P_pred_f2 = s_pred_f2.covar;
S_f2      = H_f2 * P_pred_f2 * H_f2' + R_f2;
gate_prob_f2   = 0.95;
gate_thresh_f2 = chi2inv(gate_prob_f2, 2);

% 10 random measurement states scattered around the prediction
n_meas_f2 = 10;
x_rand_f2 = mvnrnd(s_pred_f2.state', 5e4 * eye(ss_f2.num_states), n_meas_f2)';  % 4x10
z_rand_f2 = zeros(2, n_meas_f2);
is_valid_f2 = false(1, n_meas_f2);
for jj = 1:n_meas_f2
    s_tmp = tracker.makeState(ss_f2, 4*dt_f2, x_rand_f2(:,jj));
    z_rand_f2(:,jj) = msmt_f2.z_fun(s_tmp);
    nu = z_rand_f2(:,jj) - z_pred_f2;
    is_valid_f2(jj) = (nu' / S_f2 * nu) <= gate_thresh_f2;
end

% --- Figure 9.2a: x/y position space ---
fig_f2a = figure;
hold on; grid on;
% Track history
all_pos_f2 = sv_f2(ss_f2.pos_idx, :);
plot(all_pos_f2(1,:)/1e3, all_pos_f2(2,:)/1e3, 'o-', 'Color', colors(1,:), ...
     'MarkerSize', 6, 'MarkerFaceColor', colors(1,:), 'DisplayName', 'Track History');
% Predicted state + gate ellipse
xp_f2 = s_pred_f2.state(ss_f2.pos_idx);
Pp_f2 = P_pred_f2(ss_f2.pos_idx, ss_f2.pos_idx);
plot(xp_f2(1)/1e3, xp_f2(2)/1e3, 's', 'Color', colors(1,:), ...
     'MarkerSize', 6, 'MarkerFaceColor', colors(1,:), 'LineWidth', 1.5, 'DisplayName', 'Predicted State');
% Dashed line from last track state to predicted state
h_pred_line = plot([all_pos_f2(1,end), xp_f2(1)]/1e3, [all_pos_f2(2,end), xp_f2(2)]/1e3, ...
                   '--', 'Color', colors(1,:), 'LineWidth', 1.5);
utils.excludeFromLegend(h_pred_line);
theta_ell = linspace(0, 2*pi, 101);
[V_f2, D_f2] = eig(Pp_f2);
semi_f2 = sqrt(gate_thresh_f2 * max(diag(D_f2), 0)) / 1e3;
xy_ell = V_f2 * diag(semi_f2) * [cos(theta_ell); sin(theta_ell)];
h_ell = plot(xp_f2(1)/1e3 + xy_ell(1,:), xp_f2(2)/1e3 + xy_ell(2,:), '--', ...
             'Color', colors(1,:), 'DisplayName', 'Acceptance Gate');
utils.excludeFromLegend(h_ell);
% Sensors
plot(x_tdoa_f2(1,:)/1e3, x_tdoa_f2(2,:)/1e3, 'ko', 'MarkerSize', 6, ...
     'MarkerFaceColor', 'k', 'Clipping', 'off', 'DisplayName', 'TDOA Sensors');
% Random measurement source positions
valid_pos_f2 = x_rand_f2(ss_f2.pos_idx, is_valid_f2);
inval_pos_f2 = x_rand_f2(ss_f2.pos_idx, ~is_valid_f2);
if any(is_valid_f2)
    plot(valid_pos_f2(1,:)/1e3, valid_pos_f2(2,:)/1e3, 'v', 'Color', colors(3,:), ...
         'MarkerSize', 6, 'MarkerFaceColor', colors(3,:), 'DisplayName', 'Valid Measurements');
end
if any(~is_valid_f2)
    plot(inval_pos_f2(1,:)/1e3, inval_pos_f2(2,:)/1e3, '^', 'Color', colors(4,:), ...
         'MarkerSize', 6, 'MarkerFaceColor', colors(4,:), 'DisplayName', 'Invalid Measurements');
end
xlabel('x [km]'); ylabel('y [km]');
% title('Track, Prediction, and Measurements');
legend('Location','best');
utils.setPlotStyle(gca, {'widescreen'});

% --- Figure 9.2b: zeta (TDOA) space ---
fig_f2b = figure;
hold on; grid on;
% Gate ellipse in zeta space
[V_S_f2, D_S_f2] = eig(S_f2);
semi_S_f2 = sqrt(gate_thresh_f2 * max(diag(D_S_f2), 0));
xy_ell_z = V_S_f2 * diag(semi_S_f2) * [cos(theta_ell); sin(theta_ell)];
plot(z_pred_f2(1) + xy_ell_z(1,:), z_pred_f2(2) + xy_ell_z(2,:), '-', ...
     'Color', colors(1,:), 'DisplayName', 'Acceptance Gate');
plot(z_pred_f2(1), z_pred_f2(2), 'o', 'Color', colors(1,:), ...
     'MarkerSize', 6, 'MarkerFaceColor', colors(1,:), 'LineWidth', 1.5, 'DisplayName', 'Prediction');
valid_z_f2 = z_rand_f2(:,  is_valid_f2);
inval_z_f2 = z_rand_f2(:, ~is_valid_f2);
if any(is_valid_f2)
    plot(valid_z_f2(1,:), valid_z_f2(2,:), 'v', 'Color', colors(3,:), ...
         'MarkerSize', 6, 'MarkerFaceColor', colors(3,:), 'DisplayName', 'Valid Measurements');
end
if any(~is_valid_f2)
    plot(inval_z_f2(1,:), inval_z_f2(2,:), '^', 'Color', colors(4,:), ...
         'MarkerSize', 6, 'MarkerFaceColor', colors(4,:), 'DisplayName', 'Invalid Measurements');
end
xlabel('$\tau_{0,1}$ [m]'); ylabel('$\tau_{0,2}$ [m]');
% title('Prediction and Measurements in Zeta-Space');
legend('Location','best');
utils.setPlotStyle(gca, {'widescreen'});

utils.exportPlot(fig_f2a, [prefix '2a']);
utils.exportPlot(fig_f2b, [prefix '2b']);


%% Figure 9.3 - Association in Position and Measurement Space
%
% Three 2D AOA tracks (same scenario as Examples 9.1-9.3) are predicted
% forward 5 s.  Random states (3 near each predicted track + 15 background
% clutter) are generated and converted to AOA measurements.  NN association
% is run and the result is shown in both x/y space (position) and the
% two-sensor azimuth measurement space.

fprintf('Generating Figure 9.3...\n');

% Motion model and measurement model (same parameters as Examples 9.1-9.3)
mm_f3  = tracker.makeMotionModel('cv', 2, diag([25, 25]));
ss_f3  = mm_f3.state_space;
x_aoa_f3   = [750, 300;   % sensor x
              200, 800];  % sensor y
sigma_psi_f3 = 3 * pi/180;
R_f3    = sigma_psi_f3^2 * eye(2);
msmt_f3 = tracker.makeMeasurementModel(x_aoa_f3, [], [], [], [], [], R_f3);

t_msmt_f3 = 5;     % [s]
gate_prob_f3 = 0.75;

% Three pre-initialised tracks (note: track 2 uses vx=70 for this figure)
sv_f3 = {[0;    2e3;  0;   100], ...
         [1e3;  2e3;  70;  -10], ...
         [1e3;  1.3e3; 70;  50]};

P0_f3 = {[1e4,  1e2,  0,    0;
           1e2,  2e4,  0,    0;
           0,    0,    2e2,  0;
           0,    0,    0,    2e2], ...
          [4e4, -1e4,  0,    0;
          -1e4,  1e4,  0,    0;
           0,    0,    1e3, -5e2;
           0,    0,   -5e2,  5e2], ...
          [2e4,  1e4,  0,    0;
           1e4,  4e4,  0,    0;
           0,    0,    6e2,  4e2;
           0,    0,    4e2,  6e2]};

tracks_f3 = cell(3,1);
for kk = 1:3
    s = tracker.makeState(ss_f3, 0, sv_f3{kk}, P0_f3{kk});
    tracks_f3{kk} = tracker.makeTrack(s, kk);
end

% Predict to t_msmt_f3
s_pred_f3 = cell(3,1);
for kk = 1:3
    s_pred_f3{kk} = tracker.predictState(tracker.currState(tracks_f3{kk}), t_msmt_f3, mm_f3);
end

% Generate random measurement states: 3 near each predicted track + 15 background
num_near_f3 = 3;
rand_states_f3 = {};
for kk = 1:3
    mu = s_pred_f3{kk}.state';
    C  = s_pred_f3{kk}.covar;
    x_near = mvnrnd(mu, C, num_near_f3);     % num_near x 4
    for jj = 1:num_near_f3
        rand_states_f3{end+1} = tracker.makeState(ss_f3, t_msmt_f3, x_near(jj,:)'); %#ok<SAGROW>
    end
end
mu_bg_f3 = [1e3, 1e3, 0, 0];
C_bg_f3  = diag([5e4, 5e4, 100, 100]);
x_bg_f3  = mvnrnd(mu_bg_f3, C_bg_f3, 15);  % 15 x 4
for jj = 1:15
    rand_states_f3{end+1} = tracker.makeState(ss_f3, t_msmt_f3, x_bg_f3(jj,:)'); %#ok<SAGROW>
end
n_rand_f3 = numel(rand_states_f3);

% Convert random states to AOA measurements
msmts_f3 = cell(n_rand_f3, 1);
z_all_f3  = zeros(2, n_rand_f3);
for jj = 1:n_rand_f3
    z_all_f3(:,jj)  = msmt_f3.z_fun(rand_states_f3{jj});
    msmts_f3{jj}    = tracker.makeMeasurement(msmt_f3, rand_states_f3{jj}, t_msmt_f3);
end

% Predicted measurements for each track
z_pred_f3 = zeros(2, 3);
S_f3      = zeros(2, 2, 3);
for kk = 1:3
    z_pred_f3(:,kk) = msmt_f3.z_fun(s_pred_f3{kk});
    H_k = msmt_f3.h_fun(s_pred_f3{kk});
    if ndims(H_k) == 3; H_k = H_k(:,:,1)'; end
    S_f3(:,:,kk) = H_k * s_pred_f3{kk}.covar * H_k' + R_f3;
end

% NN association
[trk_idx_f3, msmt_idx_f3, ~] = tracker.associateTracks(tracks_f3, msmts_f3, ...
    t_msmt_f3, mm_f3, msmt_f3, gate_prob_f3, 'nn');

% --- Figure 9.3a: x/y space ---
fig_f3a = figure;
hold on; grid on;
scale_f3 = 1e3;
for kk = 1:3
    c = colors(kk,:);
    % Initial track position
    x0 = sv_f3{kk}(ss_f3.pos_idx);
    xp = s_pred_f3{kk}.state(ss_f3.pos_idx);
    Pp = s_pred_f3{kk}.covar(ss_f3.pos_idx, ss_f3.pos_idx);
    plot([x0(1), xp(1)]/scale_f3, [x0(2), xp(2)]/scale_f3, 'o-', 'Color', c, ...
         'DisplayName', sprintf('Track %d (pred.)', kk),'MarkerSize',6,'MarkerFaceColor',c);
    % Covariance ellipse at prediction
    [V_k, D_k] = eig(Pp);
    semi_k = sqrt(chi2inv(gate_prob_f3, 2) * max(diag(D_k), 0)) / scale_f3;
    xy_k   = V_k * diag(semi_k) * [cos(theta_ell); sin(theta_ell)];
    h_k = plot(xp(1)/scale_f3 + xy_k(1,:), xp(2)/scale_f3 + xy_k(2,:), '--', 'Color', c);
    utils.excludeFromLegend(h_k);
end
% All measurement source positions
all_rand_pos_f3 = cell2mat(cellfun(@(s) s.state(ss_f3.pos_idx), rand_states_f3, ...
                                   'UniformOutput', false));
h_msmt_f3a = plot(all_rand_pos_f3(1,:)/scale_f3, all_rand_pos_f3(2,:)/scale_f3, 'v', ...
     'Color', [0.5 0.5 0.5], 'MarkerSize', 6, 'MarkerFaceColor', [0.5 0.5 0.5], 'DisplayName', 'Measurements');
uistack(h_msmt_f3a, 'bottom');
% Sensors
plot(x_aoa_f3(1,:)/scale_f3, x_aoa_f3(2,:)/scale_f3, 'ko', ...
     'MarkerSize', 6, 'MarkerFaceColor', 'k', 'Clipping', 'off', 'DisplayName', 'DF Sensors');
xlabel('x [km]'); ylabel('y [km]');
% title('Predicted Track States with New Measurements');
legend('Location','best');
xlim([-0.5, 2.5]); ylim([-0.5, 3.5]);
utils.setPlotStyle(gca, {'widescreen'});

% --- Figure 9.3b: azimuth measurement (zeta) space ---
fig_f3b = figure;
hold on; grid on;
% All measurements as dark-gray triangles (background)
h_msmt_f3b = plot(z_all_f3(1,:), z_all_f3(2,:), 'v', 'Color', [0.5 0.5 0.5], ...
     'MarkerSize', 6, 'MarkerFaceColor', [0.5 0.5 0.5], 'DisplayName', 'Measurements');
uistack(h_msmt_f3b, 'bottom');
% For each track: predicted measurement, acceptance gate, association line, selected measurement
gate_thresh_f3 = chi2inv(gate_prob_f3, 2);
for kk = 1:3
    c  = colors(kk,:);
    zp = z_pred_f3(:,kk);
    Sk = S_f3(:,:,kk);
    % Gate ellipse
    [V_Sk, D_Sk] = eig(Sk);
    semi_Sk = sqrt(gate_thresh_f3 * max(diag(D_Sk), 0));
    xy_Sk   = V_Sk * diag(semi_Sk) * [cos(theta_ell); sin(theta_ell)];
    plot(zp(1) + xy_Sk(1,:), zp(2) + xy_Sk(2,:), '--', 'Color', c, ...
         'DisplayName', sprintf('Track %d gate', kk));
    % Predicted measurement
    plot(zp(1), zp(2), 'o', 'Color', c, 'MarkerSize', 6, 'MarkerFaceColor', c, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Track %d pred.', kk));
    % Association line + selected measurement
    mi_k = msmt_idx_f3(trk_idx_f3 == kk);
    if ~isempty(mi_k) && mi_k > 0
        za_k = msmts_f3{mi_k}.zeta;
        h_line = plot([zp(1), za_k(1)], [zp(2), za_k(2)], '-', 'Color', c);
        utils.excludeFromLegend(h_line);
        plot(za_k(1), za_k(2), 'v', 'Color', c, 'MarkerSize', 6, 'MarkerFaceColor', c, ...
             'DisplayName', sprintf('Track %d assoc.', kk));
    end
end
xlabel('$\theta_0$ [rad]'); ylabel('$\theta_1$ [rad]');
% title('Predicted Measurements and Associated Measurements');
legend('Location','best');
utils.setPlotStyle(gca, {'widescreen'});

utils.exportPlot(fig_f3a, [prefix '3a']);
utils.exportPlot(fig_f3b, [prefix '3b']);


%% Figures 9.4 and 9.5, Example 9.1 (NN Association)
fprintf('Generating Figures 9.4 and 9.5 (Example 9.1)...\n');

figs = book2_ex9_1;

utils.exportPlot(figs(1), [prefix '4']);
utils.exportPlot(figs(2), [prefix '5a']);
utils.exportPlot(figs(3), [prefix '5b']);


%% Figure 9.6, Example 9.2 (GNN Association)
fprintf('Generating Figure 9.6 (Example 9.2)...\n');

figs = book2_ex9_2;

utils.exportPlot(figs(1), [prefix '6']);


%% Figure 9.7, Example 9.3 (PDAF Association)
fprintf('Generating Figure 9.7 (Example 9.3)...\n');

figs = book2_ex9_3;

utils.exportPlot(figs(1), [prefix '7']);


%% Figure 9.8 - Covariance Growth with Missed Detections
%
% One 4-state CV track; coast for 3 consecutive missed detections and
% show how the position uncertainty ellipse grows at each step.

fprintf('Generating Figure 9.8...\n');

mm_f8  = tracker.makeMotionModel('cv', 2, 1000*eye(2));
ss_f8  = mm_f8.state_space;
dt_f8  = 1.0;    % [s] time step

sv_f8 = 1e3 * [0.00, 0.50, 1.00, 1.25;
                0.00, 0.50, 0.50, 1.00;
                1.00, 1.00, 0.50, 0.50;
                1.00, 0.00, 1.00, 0.50];
P_last_f8 = diag([2.0, 0.75, 0.4, 0.4]) * 1e4;

trk_f8 = tracker.makeTrack(tracker.makeState(ss_f8, 0*dt_f8, sv_f8(:,1)), '0');
for kk = 2:3
    trk_f8 = tracker.appendTrack(trk_f8, tracker.makeState(ss_f8, (kk-1)*dt_f8, sv_f8(:,kk)));
end
trk_f8 = tracker.appendTrack(trk_f8, tracker.makeState(ss_f8, 3*dt_f8, sv_f8(:,4), P_last_f8));

fig_f8 = figure;
hold on; grid on;

% Plot initial track history
all_pos_f8 = sv_f8(ss_f8.pos_idx, :);
plot(all_pos_f8(1,:)/1e3, all_pos_f8(2,:)/1e3, 'o-', 'Color', colors(1,:), ...
     'MarkerSize',6,'MarkerFaceColor',colors(1,:),...
     'LineWidth', 1.5, 'DisplayName', 'Initial Track');
% Covariance ellipse at final initial state (75% confidence, matching Python default)
cov_scale_f8 = sqrt(chi2inv(0.75, 2));
xc_f8 = sv_f8(ss_f8.pos_idx, end);
Pc_f8 = P_last_f8(ss_f8.pos_idx, ss_f8.pos_idx);
[V_f8, D_f8] = eig(Pc_f8);
semi_f8 = cov_scale_f8 * sqrt(max(diag(D_f8), 0)) / 1e3;
xy_ell_f8 = V_f8 * diag(semi_f8) * [cos(theta_ell); sin(theta_ell)];
h_cov = plot(xc_f8(1)/1e3 + xy_ell_f8(1,:), xc_f8(2)/1e3 + xy_ell_f8(2,:), ...
             '--', 'Color', colors(1,:), 'DisplayName', 'State Error Covariance');
utils.excludeFromLegend(h_cov);

% Coast 3 times (missed detections)
num_missed_f8 = 3;
s_prev_f8 = tracker.currState(trk_f8);
for mm_idx = 1:num_missed_f8
    c_miss = colors(mm_idx + 1, :);
    new_t_f8 = s_prev_f8.time + dt_f8;
    s_coast  = tracker.predictState(s_prev_f8, new_t_f8, mm_f8);
    trk_f8   = tracker.appendTrack(trk_f8, s_coast, true);

    % Connection line from previous to coasted position
    xp_prev = s_prev_f8.state(ss_f8.pos_idx);
    xp_new  = s_coast.state(ss_f8.pos_idx);
    h_line  = plot([xp_prev(1), xp_new(1)]/1e3, [xp_prev(2), xp_new(2)]/1e3, ...
                   '--', 'Color', c_miss);
    utils.excludeFromLegend(h_line);

    % Coasted state marker
    plot(xp_new(1)/1e3, xp_new(2)/1e3, 's', 'Color', c_miss, 'MarkerSize', 6, 'MarkerFaceColor', c_miss, ...
         'LineWidth', 1, 'DisplayName', sprintf('%d Missed Detection%s', mm_idx, ...
         char('s' * (mm_idx > 1))));

    % Growing covariance ellipse (75% confidence, matching Python default)
    Pc_miss = s_coast.covar(ss_f8.pos_idx, ss_f8.pos_idx);
    [V_m, D_m] = eig(Pc_miss);
    semi_m = cov_scale_f8 * sqrt(max(diag(D_m), 0)) / 1e3;
    xy_m   = V_m * diag(semi_m) * [cos(theta_ell); sin(theta_ell)];
    h_ell_m = plot(xp_new(1)/1e3 + xy_m(1,:), xp_new(2)/1e3 + xy_m(2,:), ...
                   '--', 'Color', c_miss);
    if mm_idx == 1
        h_ell_m.DisplayName = 'State Error Covariance';
    else
        utils.excludeFromLegend(h_ell_m);
    end

    s_prev_f8 = s_coast;
end

xlabel('x [km]'); ylabel('y [km]');
% title('Error Covariance Growth as Missed Detections Accumulate');
legend('Location','best');
utils.setPlotStyle(gca, {'widescreen'});

utils.exportPlot(fig_f8, [prefix '8']);


%% Figure 9.9, Example 9.4 (Multi-Target TDOA Tracking)
fprintf('Generating Figure 9.9 (Example 9.4)...\n');

figs = book2_ex9_4;

utils.exportPlot(figs(1), [prefix '9a']);
utils.exportPlot(figs(2), [prefix '9b']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;
