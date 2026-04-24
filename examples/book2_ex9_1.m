function figs = book2_ex9_1()
% figs = book2_ex9_1()
%
% Executes Example 9.1 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% Nearest-Neighbour (NN) Track-to-Measurement Association.
%
% Three pre-initialised 2-D AOA tracks are predicted 5 seconds forward and
% associated with four new measurements using the NN algorithm.  Three
% figures are produced:
%
%   fig1  Before association: predicted track states (with 1-sigma
%         covariance ellipses) and LS-estimated measurement positions.
%   fig2  NN assignment: dashed lines connecting matched predicted-track
%         and measurement-LS positions; hollow circles = predicted,
%         crosses = updated.
%   fig3  Final updated track states with 1-sigma covariance ellipses.
%
% Nicholas O'Donoughue
% June 2025

fprintf('Example 9.1 (NN Association)...\n');

%% ---- Scenario setup -------------------------------------------------------
[mm, msmt_model, tracks, msmts, x_aoa, R] = init_scenario();
ss        = mm.state_space;
t_msmt    = 5;      % [s] time of new measurements
gate_prob = 0.75;
scale     = 1e3;    % m -> km for plotting
colors    = get(0, 'DefaultAxesColorOrder');

%% ---- Predict all tracks to t_msmt -----------------------------------------
num_tracks = numel(tracks);
s_pred = cell(num_tracks, 1);
for kk = 1:num_tracks
    s = tracker.currState(tracks{kk});
    s_pred{kk} = tracker.predictState(s, t_msmt, mm);
end

%% ---- LS position estimates from measurements (used in figures 1 and 2) ----
num_msmts = numel(msmts);
x_ls = zeros(2, num_msmts);
P_ls = zeros(2, 2, num_msmts);
x_init_ls = [500; 1000];   % [m] initial guess for LS solver
for jj = 1:num_msmts
    x_ls(:, jj)    = triang.lsSoln(x_aoa, msmts{jj}.zeta, R, x_init_ls);
    P_ls(:, :, jj) = triang.computeCRLB(x_aoa, x_ls(:, jj), R);
end

%% ---- Figure 1: Before association -----------------------------------------
fig1 = figure;
ax1  = axes(fig1);
hold(ax1, 'on');  grid(ax1, 'on');

for kk = 1:num_tracks
    c  = colors(kk, :);
    xc = tracker.currState(tracks{kk}).state(ss.pos_idx);
    xp = s_pred{kk}.state(ss.pos_idx);
    Pp = s_pred{kk}.covar(ss.pos_idx, ss.pos_idx);
    plot(ax1, xc(1)/scale, xc(2)/scale, 'v', 'Color', c, ...
         'MarkerSize', 6, 'MarkerFaceColor', c, 'DisplayName', sprintf('Track %d (curr.)', kk));
    h_dp = plot(ax1, [xc(1), xp(1)]/scale, [xc(2), xp(2)]/scale, '--', 'Color', c, 'LineWidth', 1.5);
    utils.excludeFromLegend(h_dp);
    plot(ax1, xp(1)/scale, xp(2)/scale, 'o', 'Color', c, ...
         'MarkerSize', 6, 'MarkerFaceColor', c, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Track %d (pred.)', kk));
    h = plot_ellipse(ax1, xp, Pp, sqrt(chi2inv(gate_prob, 2)), scale);
    h.Color = c;  h.LineStyle = '--';
    utils.excludeFromLegend(h);
end

gray = [0.5 0.5 0.5];
for jj = 1:num_msmts
    hdl = plot(ax1, x_ls(1,jj)/scale, x_ls(2,jj)/scale, '^', ...
            'Color', gray, 'MarkerSize', 6, 'MarkerFaceColor', gray, 'DisplayName', 'Msmts');
    if jj>1, utils.excludeFromLegend(hdl); end
    h = plot_ellipse(ax1, x_ls(:,jj), P_ls(:,:,jj), sqrt(chi2inv(gate_prob, 2)), scale);
    h.Color = gray;  h.LineStyle = ':';
    utils.excludeFromLegend(h);
end

plot(ax1, x_aoa(1,:)/scale, x_aoa(2,:)/scale, 'ko', 'MarkerSize', 6, ...
     'MarkerFaceColor', 'k', 'Clipping', 'off', 'DisplayName', 'DF Sensors');
format_axes(ax1, 'Predicted Track States and New Measurements (Before Association)');

%% ---- NN Association -------------------------------------------------------
[trk_idx, msmt_idx, ~] = tracker.associateTracks(tracks, msmts, t_msmt, ...
                              mm, msmt_model, gate_prob, 'nn');

%% ---- Print NN distance table ----------------------------------------------
n_meas_dim = numel(msmts{1}.zeta);
num_msmts  = numel(msmts);
D = zeros(num_tracks, num_msmts);
for ti = 1:num_tracks
    for mj = 1:num_msmts
        d = tracker.computeDistance(s_pred{ti}, msmts{mj}.zeta, msmt_model);
        D(ti, mj) = d / n_meas_dim;
    end
end
sel = zeros(num_tracks, 1);
for kk = 1:numel(trk_idx)
    sel(trk_idx(kk)) = msmt_idx(kk);
end
total_cost = 0;
for kk = 1:numel(trk_idx)
    if msmt_idx(kk) > 0
        total_cost = total_cost + D(trk_idx(kk), msmt_idx(kk));
    end
end
col_w = 9;
sep   = ['+-------+', repmat([repmat('-',1,col_w),'+'], 1, num_msmts)];
fprintf('\nNN Association Distances (normalised Mahal.^2, * = selected):\n');
fprintf('%s\n', sep);
fprintf('| Track |');
for mj = 1:num_msmts
    fprintf(' %-*s|', col_w-1, sprintf('Msmt %c', 'A'+mj-1));
end
fprintf('\n%s\n', sep);
for ti = 1:num_tracks
    fprintf('| %-5d |', ti);
    for mj = 1:num_msmts
        star = '';  if sel(ti) == mj; star = '*'; end
        fprintf(' %-*s|', col_w-1, sprintf('%.2f%s', D(ti,mj), star));
    end
    fprintf('\n');
end
fprintf('%s\n', sep);
fprintf('  Total assignment cost: %.2f\n\n', total_cost);

s_upd = s_pred;   % start from predicted; overwrite assigned ones
for kk = 1:numel(trk_idx)
    ti = trk_idx(kk);
    mi = msmt_idx(kk);
    if mi > 0
        s_upd{ti} = tracker.ekfUpdate(s_pred{ti}, msmts{mi});
    end
    tracks{ti} = tracker.appendTrack(tracks{ti}, s_upd{ti}, (mi == 0));
end

%% ---- Figure 2: NN assignment visualisation --------------------------------
fig2 = figure;
ax2  = axes(fig2);
hold(ax2, 'on');  grid(ax2, 'on');

for kk = 1:num_tracks
    c  = colors(kk, :);  % color
    xc = tracks{kk}.states{end-1}.state(ss.pos_idx); % current
    xp = s_pred{kk}.state(ss.pos_idx);               % predicted
    xu = s_upd{kk}.state(ss.pos_idx);                % updated
    mi = msmt_idx(kk);
    xa = x_ls(:,mi);                                 % state est. from the assigned msmt
    Pp = s_pred{kk}.covar(ss.pos_idx, ss.pos_idx);   % predicted err
    Pu = s_upd{kk}.covar(ss.pos_idx, ss.pos_idx);    % updated err

    % Plot current track state and solid line to predicted state
    plot(ax2, xc(1)/scale, xc(2)/scale, 'v', 'Color', c, ...
         'MarkerSize', 6, 'MarkerFaceColor', c, 'DisplayName', sprintf('Track %d', kk));
    h = plot(ax2, [xc(1), xp(1)]/scale, [xc(2), xp(2)]/scale, '-', 'Color', c, 'LineWidth', 1.5);
    utils.excludeFromLegend(h);
    % Plot predicted state with error ellipse
    plot(ax2, xp(1)/scale, xp(2)/scale, 'o', 'Color', c, ...
         'MarkerSize', 6, 'MarkerFaceColor', c, 'LineWidth', 1.5, ...
         'DisplayName', 'Predicted State');
    h = plot_ellipse(ax2, xp, Pp, sqrt(chi2inv(gate_prob, 2)), scale);
    h.Color = c;  h.LineStyle = '--';
    utils.excludeFromLegend(h);
    % Dashed line from predicted state to associated msmt (no marker; it's
    % added later)
    h = plot(ax2, [xp(1), xa(1)]/scale, [xp(2), xa(2)]/scale, '--', 'Color', c, 'LineWidth', 1.5);
    utils.excludeFromLegend(h);    
end

for jj = 1:num_msmts
    hdl = plot(ax2, x_ls(1,jj)/scale, x_ls(2,jj)/scale, '^', ...
          'Color', gray, 'MarkerSize', 6, 'MarkerFaceColor', gray, 'DisplayName', 'Msmts');
    uistack(hdl, 'bottom');
    if jj>1
        utils.excludeFromLegend(hdl);
    end
end

format_axes(ax2, 'NN Association: Predicted and Updated Track States');

%% ---- Figure 3: Final updated tracks ---------------------------------------
fig3 = figure;
ax3  = axes(fig3);
hold(ax3, 'on');  grid(ax3, 'on');

for kk = 1:num_tracks
    c  = colors(kk, :);
    s  = tracker.currState(tracks{kk});
    xp = s.state(ss.pos_idx);
    Pp = s.covar(ss.pos_idx, ss.pos_idx);
    xc = tracks{kk}.states{end-1}.state(ss.pos_idx);
    h = plot(ax3, xc(1)/scale, xc(2)/scale, 'v', 'Color', c, ...
         'MarkerSize', 6, 'MarkerFaceColor', c, 'DisplayName', sprintf('Track %d', kk));
    utils.excludeFromLegend(h);
    h = plot(ax3, [xc(1), xp(1)]/scale, [xc(2), xp(2)]/scale, '-', 'Color', c, 'LineWidth', 1.5);
    utils.excludeFromLegend(h);
    plot(ax3, xp(1)/scale, xp(2)/scale, 'v', 'Color', c, ...
         'MarkerSize', 6, 'MarkerFaceColor', c, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Track %d', kk));
    h = plot_ellipse(ax3, xp, Pp, sqrt(chi2inv(gate_prob, 2)), scale);
    h.Color = c;  h.LineStyle = '-';
    utils.excludeFromLegend(h);
end

plot(ax3, x_aoa(1,:)/scale, x_aoa(2,:)/scale, 'ko', 'MarkerSize', 6, ...
     'MarkerFaceColor', 'k', 'Clipping', 'off', 'DisplayName', 'DF Sensors');
format_axes(ax3, 'Updated Trackers after NN Association');

figs = [fig1, fig2, fig3];
end


%% ---- Shared scenario initialisation --------------------------------------
function [mm, msmt_model, tracks, msmts, x_aoa, R] = init_scenario()
% AOA sensor positions [x; y] for 2 sensors (columns)
x_aoa = [750, 300;
         200, 800];

sigma_psi = 3 * pi/180;         % 3 deg AOA noise [rad]
R         = sigma_psi^2 * eye(2);

% CV motion model, 2-D, process noise sigma_a = 5 m/s^2
mm = tracker.makeMotionModel('cv', 2, diag([25, 25]));
ss = mm.state_space;

% Measurement model: pure AOA (azimuth only, 2 sensors -> 2-D zeta)
msmt_model = tracker.makeMeasurementModel(x_aoa, [], [], [], [], [], R);

% Pre-initialised track states at t = 0  [px; py; vx; vy]  (m, m/s)
state_vecs = {[0;    2e3;  0;   100], ...
              [1e3;  2e3;  80;  -10], ...
              [1e3;  1.3e3; 70;  50]};

state_covars = {[1e4,  1e2,  0,    0;
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

tracks = cell(3, 1);
for kk = 1:3
    s = tracker.makeState(ss, 0, state_vecs{kk}, state_covars{kk});
    tracks{kk} = tracker.makeTrack(s, kk);
end

% Measurements at t = 5 s  (zeta = [az_sensor1; az_sensor2], radians)
zeta_vals = {[1.811; 1.652], ...
             [1.253; 0.803], ...
             [1.140; 0.726], ...
             [1.679; 1.454]};

msmts = cell(4, 1);
for jj = 1:4
    msmts{jj} = tracker.makeMeasurement(5, zeta_vals{jj}, msmt_model);
end
end


%% ---- Plotting helpers -----------------------------------------------------
function h = plot_ellipse(ax, center_m, P_m2, n_sigma, scale)
% Plot a 1-sigma (or n_sigma) covariance ellipse.  Axis units = center_m/scale.
if nargin < 5; scale = 1; end
theta  = linspace(0, 2*pi, 101);
[V, D] = eig(P_m2);
semi   = n_sigma * sqrt(max(diag(D), 0));
xy     = V * diag(semi) * [cos(theta); sin(theta)];
h = plot(ax, (center_m(1) + xy(1,:)) / scale, ...
             (center_m(2) + xy(2,:)) / scale);
end

function format_axes(ax, title_str)
xlabel(ax, 'x [km]');
ylabel(ax, 'y [km]');
% title(ax, title_str);
legend(ax, 'Location', 'best');
xlim(ax, [-0.5, 2.0]);
ylim(ax, [-0.5, 3.0]);
utils.setPlotStyle(gca,{'widescreen'});
end
