function figs = book2_ex9_4()
% figs = book2_ex9_4()
%
% Executes Example 9.4 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% Multi-Target Tracking with False Alarms.
%
% Three aircraft fly through a 900-second scenario tracked by four ground-
% based TDOA sensors.  At each scan, the three noisy truth TDOA measurements
% are shuffled and handed to two independent trackers:
%
%   CV model  - Constant-Velocity, sigma_a = 1 m/s^2
%   CA model  - Constant-Acceleration, sigma_a = 1 m/s^2
%
% Both trackers use a TwoPointInitiator, GNN associator (gate prob = 0.95),
% MofN promoter (3-of-5), and MissedDetectionDeleter (max_missed = 3).
%
% INPUTS
%   none
%
% OUTPUTS
%   figs    [fig1, fig2, fig3]  figure handles
%           fig1  - 2x2: truth trajectories + TDOA measurement panels
%           fig2  - truth trajectories + all confirmed tracks (CV and CA)
%           fig3  - position error vs time for each target
%
% Nicholas O'Donoughue
% June 2025

fprintf('Example 9.4...\n');

%% Scenario parameters -------------------------------------------------------
t_inc    = 1;                  % s between track updates
max_time = 900;                % s total duration
ft2m     = 0.3048;
alt      = 20000 * ft2m;      % 20 000 ft -> ~6096 m

%% Target trajectories -------------------------------------------------------
[t_vec, x_tgt_1] = make_tgt_1(max_time, t_inc, alt);
[~,     x_tgt_2] = make_tgt_2(max_time, t_inc, alt);
[~,     x_tgt_3] = make_tgt_3(max_time, t_inc, alt);

x_tgts   = {x_tgt_1, x_tgt_2, x_tgt_3};
num_tgts = 3;
num_time = numel(t_vec);

%% TDOA sensor array ---------------------------------------------------------
x_tdoa  = [15e3,  0,    0, -15e3;
              0, 15e3,  0,     0;
             30,   60,  30,   60];
n_tdoa  = size(x_tdoa, 2);
ref_idx = 1;                   % first sensor as reference

sigma_toa = 1e-7;              % 100 ns timing uncertainty
C_toa     = sigma_toa^2 * eye(n_tdoa);
C_roa     = utils.constants.c^2 * C_toa;
R         = utils.resampleCovMtx(C_roa, ref_idx);
L         = chol(R, 'lower');
num_msmt  = size(R, 1);

% Least-squares solver for the TwoPointInitiator
ls_fun   = @(zeta, x0) tdoa.lsSoln(x_tdoa, zeta, C_roa, x0, [], [], [], [], ref_idx);
% CRLB function: C_roa is already in range units (m^2), so variance_is_toa=false
crlb_fun = @(x) tdoa.computeCRLB(x_tdoa, x, C_roa, ref_idx, false);

%% Motion models and measurement models --------------------------------------
sigma_a = 1;                   % m/s^2 process noise (same for both trackers)
mm_cv   = tracker.makeMotionModel('cv', 3, sigma_a^2);
mm_ca   = tracker.makeMotionModel('ca', 3, sigma_a^2);

msmt_cv = tracker.makeMeasurementModel([], x_tdoa, [], [], ref_idx, [], mm_cv.state_space, R, ls_fun, crlb_fun);
msmt_ca = tracker.makeMeasurementModel([], x_tdoa, [], [], ref_idx, [], mm_ca.state_space, R, ls_fun, crlb_fun);

%% Tracker states (independent CV and CA) ------------------------------------
ts_cv = tracker.makeTrackerState(mm_cv, msmt_cv, ...
    'gate_probability', 0.95, 'num_hits', 3, 'num_chances', 5, ...
    'max_missed', 3, 'keep_all_tracks', true);
ts_ca = tracker.makeTrackerState(mm_ca, msmt_ca, ...
    'gate_probability', 0.95, 'num_hits', 3, 'num_chances', 5, ...
    'max_missed', 3, 'keep_all_tracks', true);

%% Figure 1 setup: 2x2 panel -------------------------------------------------
scale  = 1e3;                  % plot in km
clrs   = lines(num_tgts);
num_fa = 10;                   % false alarms per scan
max_val = 30e3;                 % false-alarm extent per RDOA channel [m]

fig1   = figure;
tlo    = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax_geo = nexttile(tlo, 1);    % top-left:  truth trajectories + sensors
ax_t01 = nexttile(tlo, 2);    % top-right: TDOA sensors 0,1
ax_t02 = nexttile(tlo, 3);    % bot-left:  TDOA sensors 0,2
ax_t03 = nexttile(tlo, 4);    % bot-right: TDOA sensors 0,3
ax_msmts = [ax_t01, ax_t02, ax_t03];

% Noiseless truth trajectories (geometry panel) and TDOA curves
for ti = 1:num_tgts
    x = x_tgts{ti};
    plot(ax_geo, x(1,:)/scale, x(2,:)/scale, ...
         'Color', clrs(ti,:), 'DisplayName', sprintf('Target %d', ti));
    hold(ax_geo, 'on');

    zeta_full = tdoa.measurement(x_tdoa, x, ref_idx);   % num_msmt x num_time
    for p = 1:num_msmt
        plot(ax_msmts(p), t_vec, zeta_full(p,:)/scale, 'Color', clrs(ti,:),'DisplayName',sprintf('Target %d', ti));
        hold(ax_msmts(p), 'on');
    end
end

plot(ax_geo, x_tdoa(1,:)/scale, x_tdoa(2,:)/scale, 'ko', ...
     'MarkerFaceColor', 'k', 'MarkerSize',6,'DisplayName', 'TDOA Sensors');

% Accumulators for noisy-truth and FA scatter data added during the loop
scatter_truth_t = [];   scatter_truth_z = zeros(num_msmt, 0);
scatter_fa_t    = [];   scatter_fa_z    = zeros(num_msmt, 0);

%% Main tracking loop --------------------------------------------------------
fprintf('Running trackers across %d time steps...\n', num_time);

for idx = 1:num_time
    fprintf('.');
    t = t_vec(idx);

    % Truth measurements
    msmts = {};
    for ti = 1:num_tgts
        zeta = tdoa.measurement(x_tdoa, x_tgts{ti}(:, idx), ref_idx) + L * randn(num_msmt, 1);
        scatter_truth_t(end+1)    = t;                                  %#ok<AGROW>
        scatter_truth_z(:, end+1) = zeta;
        msmts{end+1} = tracker.makeMeasurement(t, zeta);
    end

    % False alarm measurements (uniform random in each RDOA channel)
    fa_zeta = max_val * (2*rand(num_msmt, num_fa) - 1);
    for fa_i = 1:num_fa
        msmts{end+1} = tracker.makeMeasurement(t, fa_zeta(:, fa_i));  %#ok<AGROW>
    end
    scatter_fa_t = [scatter_fa_t, repmat(t, 1, num_fa)];              %#ok<AGROW>
    scatter_fa_z = [scatter_fa_z, fa_zeta];                            %#ok<AGROW>

    % Shuffle — same permuted list handed to both trackers
    msmts = msmts(randperm(numel(msmts)));

    % Update both trackers
    ts_cv = tracker.runTrackerStep(ts_cv, msmts, t);
    ts_ca = tracker.runTrackerStep(ts_ca, msmts, t);

    fprintf(' t=%ds | CV firm=%d tent=%d | CA firm=%d tent=%d\n', ...
        round(t), numel(ts_cv.firm_tracks), numel(ts_cv.tentative_tracks), ...
        numel(ts_ca.firm_tracks), numel(ts_ca.tentative_tracks));
end

fprintf('done.\n');

% Collect all confirmed tracks (ever-promoted tracks: deleted after confirmation + still active).
% failed_tracks holds tentative tracks that never reached the promotion threshold.
all_cv = [ts_cv.deleted_tracks, ts_cv.firm_tracks];
all_ca = [ts_ca.deleted_tracks, ts_ca.firm_tracks];

fprintf('CV: %d confirmed tracks (%d deleted, %d still active, %d tentative failed)\n', ...
    numel(all_cv), numel(ts_cv.deleted_tracks), numel(ts_cv.firm_tracks), numel(ts_cv.failed_tracks));
fprintf('CA: %d confirmed tracks (%d deleted, %d still active, %d tentative failed)\n', ...
    numel(all_ca), numel(ts_ca.deleted_tracks), numel(ts_ca.firm_tracks), numel(ts_ca.failed_tracks));

%% Finish Figure 1 -----------------------------------------------------------
% Scatter the accumulated noisy-truth and FA points on the TDOA panels
msmt_titles  = {'TDOA for Sensors 0,1', 'TDOA for Sensors 0,2', 'TDOA for Sensors 0,3'};
msmt_ylabels = {'\tau_{0,1} [km]', '\tau_{0,2} [km]', '\tau_{0,3} [km]'};

for p = 1:num_msmt
    scatter(ax_msmts(p), scatter_truth_t, scatter_truth_z(p,:)/scale, ...
            3, [0 0 1], 'filled', 'DisplayName', 'Truth Msmts');
    scatter(ax_msmts(p), scatter_fa_t, scatter_fa_z(p,:)/scale, ...
            3, [0.7 0.7 0.7], 'filled', 'DisplayName', 'False Alarms');
    title(ax_msmts(p),  msmt_titles{p},  'FontSize', 10);
    xlabel(ax_msmts(p), 'Time [s]',      'FontSize', 8);
    ylabel(ax_msmts(p), msmt_ylabels{p}, 'FontSize', 8);
    legend(ax_msmts(p), 'FontSize', 8,'Location','SouthEast');
    grid(ax_msmts(p), 'on');
    ax_msmts(p).FontSize = 8;
end

title(ax_geo,  'Target Trajectories', 'FontSize', 10);
xlabel(ax_geo, 'East [km]',           'FontSize', 8);
ylabel(ax_geo, 'North [km]',          'FontSize', 8);
xlim(ax_geo, [-75, 125]);
ylim(ax_geo, [-50, 200]);
legend(ax_geo, 'FontSize', 8);
grid(ax_geo, 'on');
ax_geo.FontSize = 8;

%% Figure 2: truth trajectories + confirmed tracks ---------------------------
fig2 = figure;
ax2  = axes('Parent', fig2);
hold(ax2, 'on');

for ti = 1:num_tgts
    x = x_tgts{ti};
    plot(ax2, x(1,:)/scale, x(2,:)/scale, ...
         'Color', clrs(ti,:), 'DisplayName', sprintf('Truth %d', ti));
end

lbl = 'CV Tracks';
first_trk = true;
for trk = all_cv
    pos_idx = mm_cv.state_space.pos_idx;
    states  = trk{1}.states;
    x_trk   = cell2mat(cellfun(@(s) s.state(pos_idx), states, 'UniformOutput', false));
    hdl = plot(ax2, x_trk(1,:)/scale, x_trk(2,:)/scale, 'r--', 'DisplayName', lbl);
    if ~first_trk
        utils.excludeFromLegend(hdl);
    end
    first_trk = false;
end

lbl = 'CA Tracks';
first_trk = true;
for trk = all_ca
    pos_idx = mm_ca.state_space.pos_idx;
    states  = trk{1}.states;
    x_trk   = cell2mat(cellfun(@(s) s.state(pos_idx), states, 'UniformOutput', false));
    hdl = plot(ax2, x_trk(1,:)/scale, x_trk(2,:)/scale, 'm-.', 'DisplayName', lbl);
    if ~first_trk
        utils.excludeFromLegend(hdl);
    end
    first_trk = false;
end

scatter(ax2, x_tdoa(1,:)/scale, x_tdoa(2,:)/scale, 60, 'ks', 'filled', ...
        'DisplayName', 'TDOA Sensors');
xlim(ax2, [-75, 125]);
ylim(ax2, [-50, 200]);
xlabel(ax2, 'East [km]');
ylabel(ax2, 'North [km]');
title(ax2, 'Example 9.4: Truth Trajectories and Tracker Output (CV vs CA)');
legend(ax2, 'FontSize', 8);
grid(ax2, 'on');
utils.setPlotStyle(ax2, {'widescreen'});

%% Figure 3: position error vs time -----------------------------------------
fig3 = figure;
ax3  = axes('Parent', fig3);
hold(ax3, 'on');

for ti = 1:num_tgts
    cv_err = compute_position_error(all_cv, t_vec, x_tgts{ti}, mm_cv.state_space.pos_idx);
    ca_err = compute_position_error(all_ca, t_vec, x_tgts{ti}, mm_ca.state_space.pos_idx);
    plot(ax3, t_vec, cv_err/scale, '--', 'Color', clrs(ti,:), ...
         'DisplayName', sprintf('Target %d CV', ti));
    plot(ax3, t_vec, ca_err/scale, '-.',  'Color', clrs(ti,:), ...
         'DisplayName', sprintf('Target %d CA', ti));
end

title(ax3,  'Example 9.4: Horizontal Position Error vs Time (CV dashed, CA dash-dot)', 'FontSize', 10);
xlabel(ax3, 'Time [s]',  'FontSize', 8);
ylabel(ax3, 'Horizontal Error [km]','FontSize', 8);
legend(ax3, 'FontSize', 8);
grid(ax3, 'on');
ax3.FontSize = 8;
utils.setPlotStyle(ax3, {'widescreen'});

figs = [fig1, fig2, fig3];

end  % book2_ex9_4


%% ==========================================================================
function [t_vec, x_tgt] = make_tgt_1(max_time, t_inc, alt)
% Target 1: due east for 3 min, quarter-circle right turn (to south), then
% due south.
x0     = [-50e3; 50e3; alt];
vel    = 200;            % m/s
t_e    = 3*60;           % east-leg duration [s]
r_turn = 50e3;           % turn radius [m]
t_turn = (pi/2) * r_turn / vel;   % quarter-circle arc time [s]

t_e_vec    = (0     : t_inc : t_e)';
t_turn_vec = (t_inc : t_inc : t_turn)';
t_s_vec    = (t_inc : t_inc : max_time)';

% East leg
x_e = x0 + [vel; 0; 0] * t_e_vec';

% Turn: parameterise by arc fraction (uniform time steps)
ang    = (pi/2) * t_turn_vec / t_turn;
x_turn = x_e(:, end) + r_turn * [sin(ang)'; cos(ang)' - 1; zeros(1, numel(ang))];

% South leg
x_s = x_turn(:, end) + [0; -vel; 0] * t_s_vec';

x_full = [x_e, x_turn, x_s];
t_full = (0 : size(x_full, 2)-1)' * t_inc;

mask  = t_full <= max_time;
t_vec = t_full(mask);
x_tgt = x_full(:, mask);
end


%% ==========================================================================
function [t_vec, x_tgt] = make_tgt_2(max_time, t_inc, alt)
% Target 2: due north-east for 3 min, 180-degree left turn, then due
% south-west.
x0    = [-50e3; 75e3; alt];
vel   = 210;             % m/s
t_ne  = 3*60;            % NE-leg duration [s]
r_turn = 50e3;           % turn radius [m]
t_turn = pi * r_turn / vel;        % semicircle arc time [s]
n     = 1/sqrt(2);       % NE / SW unit-vector component

t_ne_vec   = (0     : t_inc : t_ne)';
t_turn_vec = (t_inc : t_inc : t_turn)';
t_sw_vec   = (t_inc : t_inc : max_time)';

% North-east leg
x_ne = x0 + [vel*n; vel*n; 0] * t_ne_vec';

% 180-degree turn: circle parameterised from angle 3*pi/4 to -pi/4
% (clockwise sweep, matching the NE -> SW reversal)
num_turn_pts = numel(t_turn_vec);
ang      = 3*pi/4 - pi * (0:num_turn_pts-1)' / num_turn_pts;
x_turn0  = [cos(ang)'; sin(ang)'; zeros(1, num_turn_pts)];
x_turn   = x_ne(:, end) + r_turn * (x_turn0 - x_turn0(:, 1));

% South-west leg
x_sw = x_turn(:, end) + [-vel*n; -vel*n; 0] * t_sw_vec';

x_full = [x_ne, x_turn, x_sw];
t_full = (0 : size(x_full, 2)-1)' * t_inc;

mask  = t_full <= max_time;
t_vec = t_full(mask);
x_tgt = x_full(:, mask);
end


%% ==========================================================================
function [t_vec, x_tgt] = make_tgt_3(max_time, t_inc, alt)
% Target 3: straight south-east leg for the full scenario duration.
x0  = [-50e3; 125e3; alt];
vel = 170;               % m/s
n   = 1/sqrt(2);         % SE unit-vector component

t_vec = (0 : t_inc : max_time)';
x_tgt = x0 + [vel*n; -vel*n; 0] * t_vec';
end


%% ==========================================================================
function err = compute_position_error(all_tracks, t_vec, x_truth, pos_idx)
% For each element of t_vec, return the minimum horizontal (x-y) distance
% from x_truth(1:2, idx) to the x-y position stored in any confirmed track
% at that time.  Z is excluded because the sensor array is nearly coplanar
% (z~30-60 m) while targets fly at ~6096 m, making altitude nearly
% unobservable via TDOA alone (CRLB_z >> CRLB_xy for this geometry).
% Returns NaN where no track has a state at that timestamp.

num_time = numel(t_vec);
err      = nan(1, num_time);

for t_idx = 1:num_time
    t      = t_vec(t_idx);
    x_true = x_truth(1:2, t_idx);   % horizontal truth position only
    min_d  = inf;

    for trk_j = 1:numel(all_tracks)
        states = all_tracks{trk_j}.states;
        for s_k = 1:numel(states)
            if abs(states{s_k}.time - t) < 0.5    % within half a step
                x_est = states{s_k}.state(pos_idx(1:2));   % xy only
                d = norm(x_true - x_est);
                if d < min_d
                    min_d = d;
                end
            end
        end
    end

    if isfinite(min_d)
        err(t_idx) = min_d;
    end
end
end
