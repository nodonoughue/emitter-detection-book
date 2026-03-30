function figs = book2_ex9_2()
% figs = book2_ex9_2()
%
% Executes Example 9.2 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% Global Nearest-Neighbour (GNN) Track-to-Measurement Association.
%
% Same three-track, four-measurement AOA scenario as Example 9.1, but
% association is solved as a global assignment problem (Hungarian/Munkres
% algorithm) rather than the greedy sequential NN search.  One figure is
% produced showing the final updated track states.
%
% Nicholas O'Donoughue
% June 2025

fprintf('Example 9.2 (GNN Association)...\n');

%% ---- Scenario setup -------------------------------------------------------
[mm, msmt_model, tracks, msmts, x_aoa, ~] = init_scenario();
ss        = mm.state_space;
t_msmt    = 5;      % [s]
gate_prob = 0.75;
scale     = 1e3;    % m -> km
colors    = get(0, 'DefaultAxesColorOrder');

%% ---- GNN Association ------------------------------------------------------
[trk_idx, msmt_idx, ~] = tracker.associateTracks(tracks, msmts, t_msmt, ...
                              mm, msmt_model, gate_prob, 'gnn');

for kk = 1:numel(trk_idx)
    ti = trk_idx(kk);
    mi = msmt_idx(kk);
    s      = tracker.currState(tracks{ti});
    s_pred = tracker.predictState(s, t_msmt, mm);
    if mi > 0
        s_upd = tracker.ekfUpdateState(s_pred, msmts{mi}.zeta, msmt_model);
        tracks{ti} = tracker.appendTrack(tracks{ti}, s_upd, false);
    else
        tracks{ti} = tracker.appendTrack(tracks{ti}, s_pred, true);
    end
end

%% ---- Figure: Updated tracks -----------------------------------------------
fig1 = figure;
ax   = axes(fig1);
hold(ax, 'on');  grid(ax, 'on');

for kk = 1:numel(tracks)
    c  = colors(kk, :);
    s  = tracker.currState(tracks{kk});
    xp = s.state(ss.pos_idx);
    Pp = s.covar(ss.pos_idx, ss.pos_idx);
    plot(ax, xp(1)/scale, xp(2)/scale, 'v', 'Color', c, ...
         'MarkerSize', 10, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Track %d', kk));
    h = plot_ellipse(ax, xp, Pp, 1, scale);
    h.Color = c;
    utils.excludeFromLegend(h);
end

plot(ax, x_aoa(1,:)/scale, x_aoa(2,:)/scale, 'ks', 'MarkerSize', 10, ...
     'MarkerFaceColor', 'k', 'DisplayName', 'DF Sensors');
format_axes(ax, 'Updated Trackers after GNN Association');

figs = fig1;
end


%% ---- Shared scenario initialisation --------------------------------------
function [mm, msmt_model, tracks, msmts, x_aoa, R] = init_scenario()
x_aoa = [750, 300;
         200, 800];

sigma_psi = 3 * pi/180;
R         = sigma_psi^2 * eye(2);

mm = tracker.makeMotionModel('cv', 2, diag([25, 25]));
ss = mm.state_space;

[z_fun, h_fun] = tracker.makeMeasurementModel(x_aoa, [], [], [], [], [], ss);
msmt_model = tracker.makeMsmtModel(z_fun, h_fun, R, ss);

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

zeta_vals = {[1.811; 1.652], ...
             [1.253; 0.803], ...
             [1.140; 0.726], ...
             [1.679; 1.454]};

msmts = cell(4, 1);
for jj = 1:4
    msmts{jj} = tracker.makeMeasurement(5, zeta_vals{jj});
end
end


%% ---- Plotting helpers -----------------------------------------------------
function h = plot_ellipse(ax, center_m, P_m2, n_sigma, scale)
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
title(ax, title_str);
legend(ax, 'Location', 'best');
xlim(ax, [-0.5, 2.0]);
ylim(ax, [-0.5, 3.0]);
end
