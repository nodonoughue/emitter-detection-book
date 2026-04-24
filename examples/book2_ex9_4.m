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
t_inc    = 10;                 % s between track updates
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
x_tdoa  = [ 0,  10.6e3,    0, -10.6e3;
            0, -10.6e3, 15e3, -10.6e3;
           60,      30,   30,      30];
n_tdoa  = size(x_tdoa, 2);
ref_idx = 1;                   % first sensor as reference

sigma_toa = 5e-8;              % 50 ns timing uncertainty
C_toa     = sigma_toa^2 * eye(n_tdoa);
C_roa     = utils.constants.c^2 * C_toa;
R         = utils.resampleCovMtx(C_roa, ref_idx);
L         = chol(R, 'lower');
num_msmt  = size(R, 1);

% Least-squares solver for the TwoPointInitiator
% Restrict the LS solver to 25 steps; we don't need exact positions, just
% coarse ones for initializing tracks.
max_steps = 25;
min_alt = 300; % it's flying
max_alt = 40000; % it's not in space
bnd = utils.constraints.boundedAlt(min_alt, max_alt,'flat');
ls_fun   = @(zeta, x0) tdoa.lsSolnBounded(x_tdoa, zeta, C_roa, x0, bnd, ...
    [], max_steps, [], [], ref_idx);
% CRLB function: C_roa is already in range units (m^2), so variance_is_toa=false
crlb_fun = @(x) tdoa.computeCRLB(x_tdoa, x, C_roa, ref_idx, false);

%% Motion models and measurement models --------------------------------------
q_a_cv = 6;              % m/s^2 process noise for CV tracker
q_a_ca = .3;              % m/s^2 process noise for CA tracker
mm_cv   = tracker.makeMotionModel('cv', 3, q_a_cv^2);
mm_ca   = tracker.makeMotionModel('ca', 3, q_a_ca^2);

msmt = tracker.makeMeasurementModel([], x_tdoa, [], [], ref_idx, [], R, ls_fun, crlb_fun);

%% Tracker states (independent CV and CA) ------------------------------------
target_max_vel   = 350;   % m/s  — conservative upper bound for subsonic aircraft
target_max_accel = 10;    % m/s² — generous bound for a maneuvering aircraft

% gate_probability: Python uses 0.95 for CV and CA.
% At 100+ km range the TDOA Jacobian magnitude is ~0.07, so a 2 km EKF
% position error maps to ~140 m TDOA prediction mismatch → d²≈30.
% CA needs a wide gate to avoid coasting through those steps; CV is tighter
% because its process noise (q_a=3) keeps P from growing as large.
ts_cv = tracker.makeTrackerState(mm_cv, msmt, ...
    'gate_probability', 0.95, 'num_hits', 3, 'num_chances', 5, ...
    'max_missed', 3, 'keep_all_tracks', true, ...
    'target_max_velocity', target_max_vel);
ts_ca = tracker.makeTrackerState(mm_ca, msmt, ...
    'gate_probability', 0.95, 'num_hits', 3, 'num_chances', 5, ...
    'max_missed', 3, 'keep_all_tracks', true, ...
    'target_max_velocity', target_max_vel, ...
    'target_max_acceleration', target_max_accel);

%% Figure 1 setup: 2x2 panel -------------------------------------------------
scale  = 1e3;                  % plot in km
clrs   = lines(num_tgts);
num_fa = 10;                   % false alarms per scan
max_val = 15e3;                % false-alarm extent per RDOA channel [m]

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
    hdl=plot(ax_geo, x(1,:)/scale, x(2,:)/scale, ...
             'Color', clrs(ti,:), 'DisplayName', 'Targets');
    if ti > 1, utils.excludeFromLegend(hdl); end
    hold(ax_geo, 'on');
    mid    = floor(num_time / 2);
    tip    = x(1:2, mid+1) / scale;               % arrowhead tip [km]
    base_c = x(1:2, mid)   / scale;               % arrowhead base centre [km]
    shaft  = 5 * (tip - base_c);                   % direction vector [km], scaled 5x
    perp   = [-shaft(2); shaft(1)] / norm(shaft);  % unit perpendicular
    hw     = 0.6 * norm(shaft);                    % arrowhead half-width [km]
    patch(ax_geo, ...
          [base_c(1)+shaft(1), base_c(1)+hw*perp(1), base_c(1)-hw*perp(1)], ...
          [base_c(2)+shaft(2), base_c(2)+hw*perp(2), base_c(2)-hw*perp(2)], ...
          clrs(ti,:), 'EdgeColor', 'none', 'HandleVisibility', 'off');

    zeta_full = tdoa.measurement(x_tdoa, x, ref_idx);   % num_msmt x num_time
    for p = 1:num_msmt
        hdl = plot(ax_msmts(p), t_vec, zeta_full(p,:)/scale, 'Color', clrs(ti,:),'DisplayName','Targets');
        if ti > 1, utils.excludeFromLegend(hdl); end
        hold(ax_msmts(p), 'on');
    end
end

plot(ax_geo, x_tdoa(1,:)/scale, x_tdoa(2,:)/scale, 'ko', ...
     'MarkerFaceColor', 'k', 'MarkerSize', 3, 'Clipping', 'off', 'DisplayName', 'TDOA Sensors');

% Accumulators for noisy-truth and FA scatter data (pre-allocated)
scatter_truth_t = zeros(1,        num_time * num_tgts);
scatter_truth_z = zeros(num_msmt, num_time * num_tgts);
scatter_fa_t    = zeros(1,        num_time * num_fa);
scatter_fa_z    = zeros(num_msmt, num_time * num_fa);
truth_fill = 0;
fa_fill    = 0;

%% Main tracking loop --------------------------------------------------------
fprintf('Running trackers across %d time steps...\n', num_time);

for idx = 1:num_time
    t = t_vec(idx);

    % Truth measurements
    msmts = cell(1, num_tgts + num_fa);
    for ti = 1:num_tgts
        zeta = tdoa.measurement(x_tdoa, x_tgts{ti}(:, idx), ref_idx) + L * randn(num_msmt, 1);
        truth_fill = truth_fill + 1;
        scatter_truth_t(truth_fill)      = t;
        scatter_truth_z(:, truth_fill)   = zeta;
        msmts{ti} = tracker.makeMeasurement(t, zeta, msmt);
    end

    % False alarm measurements (uniform random in each RDOA channel)
    fa_zeta = max_val * (2*rand(num_msmt, num_fa) - 1);
    fa_cols = fa_fill + (1:num_fa);
    scatter_fa_t(fa_cols)    = t;
    scatter_fa_z(:, fa_cols) = fa_zeta;
    fa_fill = fa_fill + num_fa;
    for fa_i = 1:num_fa
        msmts{num_tgts + fa_i} = tracker.makeMeasurement(t, fa_zeta(:, fa_i), msmt);
    end

    % Shuffle — same permuted list handed to both trackers
    msmts = msmts(randperm(numel(msmts)));

    % Update both trackers (each with its own measurement set / model)
    ts_cv = tracker.runTrackerStep(ts_cv, msmts, t);
    ts_ca = tracker.runTrackerStep(ts_ca, msmts, t);

    fprintf(' t=%ds | CV firm=%d tent=%d | CA firm=%d tent=%d\n', ...
        round(t), numel(ts_cv.firm_tracks), numel(ts_cv.tentative_tracks), ...
        numel(ts_ca.firm_tracks), numel(ts_ca.tentative_tracks));
end

fprintf('done.\n');
% After the main loop, before collecting all_cv / all_ca
fprintf('\n--- Final CV firm track positions [km] ---\n');
for ii = 1:numel(ts_cv.firm_tracks)
    s = tracker.currState(ts_cv.firm_tracks{ii});
    p = s.state(mm_cv.state_space.pos_idx);
    fprintf('  Track %d: [%.1f, %.1f, %.1f] km\n', ii, p(1)/1e3, p(2)/1e3, p(3)/1e3);
end
fprintf('Truth positions at t=900s [km]:\n');
for ti = 1:num_tgts
    p = x_tgts{ti}(:,end);
    fprintf('  Target %d: [%.1f, %.1f, %.1f] km\n', ti, p(1)/1e3, p(2)/1e3, p(3)/1e3);
end
% Collect all confirmed tracks (ever-promoted tracks: deleted after confirmation + still active).
% failed_tracks holds tentative tracks that never reached the promotion threshold.
all_cv = [ts_cv.deleted_tracks, ts_cv.firm_tracks];
all_ca = [ts_ca.deleted_tracks, ts_ca.firm_tracks];

fprintf('CV: %d confirmed tracks (%d deleted, %d still active, %d tentative failed)\n', ...
    numel(all_cv), numel(ts_cv.deleted_tracks), numel(ts_cv.firm_tracks), numel(ts_cv.failed_tracks));
fprintf('CA: %d confirmed tracks (%d deleted, %d still active, %d tentative failed)\n', ...
    numel(all_ca), numel(ts_ca.deleted_tracks), numel(ts_ca.firm_tracks), numel(ts_ca.failed_tracks));

%% Finish Figure 1 -----------------------------------------------------------
% Trim pre-allocated scatter arrays to actual fill size
scatter_truth_t = scatter_truth_t(1:truth_fill);
scatter_truth_z = scatter_truth_z(:, 1:truth_fill);
scatter_fa_t    = scatter_fa_t(1:fa_fill);
scatter_fa_z    = scatter_fa_z(:, 1:fa_fill);

% Scatter the accumulated noisy-truth and FA points on the TDOA panels
msmt_titles  = {'TDOA for Sensors 0,1', 'TDOA for Sensors 0,2', 'TDOA for Sensors 0,3'};
msmt_ylabels = {'$R_{0,1}$ [km]', '$R_{0,2}$ [km]', '$R_{0,3}$ [km]'};

for p = 1:num_msmt
    scatter(ax_msmts(p), scatter_truth_t, scatter_truth_z(p,:)/scale, ...
            3, [0 0 1], 'filled', 'DisplayName', 'Truth Msmts');
    scatter(ax_msmts(p), scatter_fa_t, scatter_fa_z(p,:)/scale, ...
            3, [0.7 0.7 0.7], 'filled', 'DisplayName', 'False Alarms');
    title(ax_msmts(p),  msmt_titles{p},  'FontSize', 8);
    xlabel(ax_msmts(p), 'Time [s]',      'FontSize', 6);
    ylabel(ax_msmts(p), msmt_ylabels{p}, 'FontSize', 6);
    legend(ax_msmts(p), 'FontSize', 6,'Location','SouthEast');
    grid(ax_msmts(p), 'on');
    ax_msmts(p).FontSize = 6;
end

title(ax_geo,  'Target Trajectories', 'FontSize', 8);
xlabel(ax_geo, 'East [km]',           'FontSize', 6);
ylabel(ax_geo, 'North [km]',          'FontSize', 6);
xlim(ax_geo, [-75, 125]);
ylim(ax_geo, [-50, 200]);
legend(ax_geo, 'FontSize', 6,'Location','northwest');
grid(ax_geo, 'on');
ax_geo.FontSize = 6;

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

scatter(ax2, x_tdoa(1,:)/scale, x_tdoa(2,:)/scale, 60, 'k', 'filled', ...
        'Clipping', 'off', 'DisplayName', 'TDOA Sensors');
xlim(ax2, [-75, 125]);
ylim(ax2, [-50, 200]);
xlabel(ax2, 'East [km]');
ylabel(ax2, 'North [km]');
% title(ax2, 'Example 9.4: Truth Trajectories and Tracker Output (CV vs CA)');
legend(ax2, 'FontSize', 8);
grid(ax2, 'on');
utils.setPlotStyle(ax2, {'widescreen'});

figs = [fig1, fig2];

end  % book2_ex9_4


%% ==========================================================================
function [t_vec, x_tgt] = make_tgt_1(max_time, t_inc, alt)
% Target 1: due east for 3 min, quarter-circle right turn (to south), then
% due south.
x0     = [-50e3; 80e3; alt];
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

