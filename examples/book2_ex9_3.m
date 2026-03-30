function figs = book2_ex9_3()
% figs = book2_ex9_3()
%
% Executes Example 9.3 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% Probabilistic Data Association Filter (PDAF) Track-to-Measurement
% Association.
%
% Same three-track, four-measurement AOA scenario as Examples 9.1 and 9.2,
% but association uses soft weighting (PDAF) rather than hard assignment.
% For each track, all gated measurements contribute to the EKF update via
% a weighted Gaussian mixture:
%
%   beta_j  = L_j / (p_miss + sum_k L_k)     for each gated measurement j
%   beta_0  = p_miss / (p_miss + sum_k L_k)  null hypothesis (no detection)
%   p_miss  = 1 - P_d * P_g                  (P_d = 1, P_g = gate_prob)
%
% where L_j = N(nu_j; 0, S) is the Gaussian likelihood of innovation nu_j.
% The fused state is the weighted mean of the per-hypothesis EKF updates
% (including a coasted null hypothesis), and the fused covariance is the
% corresponding mixture covariance (sum of per-hypothesis covars plus the
% spread-of-means term).
%
% One figure is produced showing the final updated track states.
%
% Nicholas O'Donoughue
% June 2025

fprintf('Example 9.3 (PDAF Association)...\n');

%% ---- Scenario setup -------------------------------------------------------
[mm, msmt_model, tracks, msmts, x_aoa, ~] = init_scenario();
ss        = mm.state_space;
t_msmt    = 5;      % [s]
gate_prob = 0.75;
scale     = 1e3;    % m -> km
colors    = get(0, 'DefaultAxesColorOrder');

%% ---- PDAF update for each track -------------------------------------------
num_tracks = numel(tracks);
for kk = 1:num_tracks
    s      = tracker.currState(tracks{kk});
    s_pred = tracker.predictState(s, t_msmt, mm);
    s_upd  = pda_update_track(s_pred, msmts, msmt_model, gate_prob);
    tracks{kk} = tracker.appendTrack(tracks{kk}, s_upd, false);
end

%% ---- Figure: Updated tracks -----------------------------------------------
fig1 = figure;
ax   = axes(fig1);
hold(ax, 'on');  grid(ax, 'on');

for kk = 1:num_tracks
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
format_axes(ax, 'Updated Trackers after PDAF Association');

figs = fig1;
end


%% ---- PDAF update ----------------------------------------------------------
function s_upd = pda_update_track(s_pred, msmts, msmt_model, gate_prob, pd)
% pda_update_track  Probabilistic Data Association EKF update for one track.
%
% For each gated measurement the standard EKF update is computed; a null
% (coast) hypothesis is always included.  The resulting state is the
% weighted Gaussian mixture of all hypotheses.
%
% INPUTS
%   s_pred      Predicted State struct
%   msmts       Cell array of all Measurement structs for the current scan
%   msmt_model  Measurement model struct from makeMsmtModel
%   gate_prob   Chi-square gate probability (e.g. 0.75)
%   pd          Detection probability (default: 1.0)

if nargin < 5; pd = 1.0; end

ss         = s_pred.state_space;
n_st       = ss.num_states;
n_all_m    = numel(msmts);
n_msmt_dim = numel(msmts{1}.zeta);

% Chi-squared gate threshold
gate_thresh = chi2inv(gate_prob, n_msmt_dim);

% Pre-compute EKF quantities at the predicted state
z_hat = msmt_model.z_fun(s_pred);
H     = msmt_model.h_fun(s_pred);
if ndims(H) == 3                          % handle 3-D Jacobian (num_src > 1)
    H = H(:, :, 1)';
end
P = s_pred.covar;
R = msmt_model.R;
S = H * P * H' + R;
K = P * H' / S;

% Find gated measurements and their Gaussian likelihoods
gated_idx  = [];
gated_like = [];
for jj = 1:n_all_m
    nu  = msmts{jj}.zeta(:) - z_hat(:);
    d2  = nu' / S * nu;          % Mahalanobis distance squared
    if d2 <= gate_thresh
        n_dim  = numel(nu);
        L_j    = exp(-0.5 * d2) / sqrt((2*pi)^n_dim * det(S));
        gated_idx(end+1)  = jj;   %#ok<AGROW>
        gated_like(end+1) = L_j;  %#ok<AGROW>
    end
end

% Null hypothesis likelihood: p_miss = 1 - Pd * Pg
p_miss = 1 - pd * gate_prob;

% Build and normalise weights  [gated_1 ... gated_M  null]
all_weights = [gated_like, p_miss];
all_weights = all_weights / sum(all_weights);

% Per-hypothesis state updates
n_gated = numel(gated_idx);
n_hyp   = n_gated + 1;          % gated measurements + null
all_states = zeros(n_st, n_hyp);
all_covars = zeros(n_st, n_st, n_hyp);

I_nst = eye(n_st);
for kk = 1:n_gated
    jj   = gated_idx(kk);
    nu_j = msmts{jj}.zeta(:) - z_hat(:);
    all_states(:, kk)    = s_pred.state + K * nu_j;
    all_covars(:, :, kk) = (I_nst - K * H) * P;
end

% Null hypothesis: coast (no measurement update)
all_states(:, end)    = s_pred.state;
all_covars(:, :, end) = P;

% Weighted mean
x_fused = all_states * all_weights(:);

% Weighted covariance + spread-of-means
P_fused = zeros(n_st);
for kk = 1:n_hyp
    delta   = all_states(:, kk) - x_fused;
    P_fused = P_fused + all_weights(kk) * (all_covars(:, :, kk) + delta * delta');
end

s_upd = tracker.makeState(ss, s_pred.time, x_fused, P_fused);
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
