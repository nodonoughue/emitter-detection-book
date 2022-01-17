% Draw Figures - Chapter 7
%
% This script generates all of the figures that appear in
% Chapter 7 of the textbook.
%
% Nicholas O'Donoughue
% 11 January 2022

% Clear Figures
close all;

% Set up directory and filename for figures
dirNm = fullfile(pwd,'figures','practical_geo');
if ~exist(dirNm,'dir')
    mkdir(dirNm);
end
prefix = fullfile(dirNm,'fig7_');

% Initialize Plot Preference
doFullColor=true;
utils.initPlotSettings;
colors=get(groot,'DefaultAxesColorOrder');

% Reset the random number generator, to ensure reproducability
rng('default') ; 

addpath('examples');

%% Figure 7.1 and 7.2, Example 7.1
figs = book2_ex7_1;

utils.exportPlot(figs(1), [prefix '1']);
utils.exportPlot(figs(2), [prefix '2a']);
utils.exportPlot(figs(3), [prefix '2b']);

%% Figure 7.3 and 7.4, Example 7.2
figs = book2_ex7_2;

utils.exportPlot(figs(1), [prefix '3']);
utils.exportPlot(figs(2), [prefix '4a']);
utils.exportPlot(figs(3), [prefix '4b']);

%% Figure 7.5, Example 7.3
figs = book2_ex7_3;

utils.exportPlot(figs(1), [prefix '5a']);
utils.exportPlot(figs(2), [prefix '5b']);

%% Figure 7.6
% Illustration of TDOA changes over time

x_tgt = [0; 0];
x_tdoa = [-5, 5, 15;
          -50, -40, -50]*1e3;
v_tdoa = [0, 0, 0;
          1e3, 1e3, 1e3];

t = 0:.1:100;

num_sensors = size(x_tdoa,2);
num_t = numel(t);

x_tdoa_full = x_tdoa + v_tdoa .* reshape(t,1,1,num_t);

% Plot Geometry
fig6a = figure;
colors = get(gca,'ColorOrder');
plot(x_tgt(1), x_tgt(2),'^', 'DisplayName','Target');
hold on;
hdl=plot(squeeze(x_tdoa_full(1,:,:))', squeeze(x_tdoa_full(2,:,:))','Color',colors(2,:),'DisplayName','TDOA Sensors');
utils.excludeFromLegend(hdl(2:end));
hdl=scatter(x_tdoa(1,:),x_tdoa(2,:),'o','filled','MarkerEdgeColor',.8*colors(2,:),'MarkerFaceColor',colors(2,:));
utils.excludeFromLegend(hdl);
legend('Location','NorthWest');

utils.setPlotStyle(gca,{'widescreen','equal'});

% Compute TDOA as a function of time
zeta = zeros(num_sensors-1,num_t);
for idx=1:num_t
    this_x = squeeze(x_tdoa_full(:,:,idx));
    zeta(:,idx) = tdoa.measurement(this_x,x_tgt,1);
end

fig6b=figure;
plot(t,zeta);
legend('TDOA_{1,2}','TDOA_{1,3}');
grid on;
xlabel('Time [s]');
ylabel('Range Difference Measurement [m]');

utils.setPlotStyle(gca,{'widescreen'});

utils.exportPlot(fig6a,[prefix '6a']);
utils.exportPlot(fig6b,[prefix '6b']);

%% Figure 7.7
figs = book2_ex7_4;

utils.exportPlot(figs(1), [prefix '7a']);
utils.exportPlot(figs(2), [prefix '7b']);

%% Figure 7.8 and 7.9 Degenerate Geometry
x_tgt = [0; 10e3];

x_init = [0; 0];
v_init = [0; 50];
a_init = [1; 0]; % m/s^2
t_turn = 100;    % s, at this point, flip the acceleration
t_full = 2*t_turn;
dt = 1;
t_vec = dt:dt:t_full;
a_aoa = cat(2,[0;0],a_init * ((t_vec < t_turn) - (t_vec > t_turn)));
v_aoa = v_init + cumsum(a_aoa * dt,2);
x_aoa = x_init + cumsum(v_aoa * dt,2);

theta_unc = 5; % +/- 5 degree uncertainty interval

fig8a=figure;
hdl_traj = plot(x_aoa(1,:), x_aoa(2,:),'DisplayName','Sensor Trajectory');
hold on;

% Draw bearings at time markers
idx_set = [1, fix(numel(t_vec)/2), numel(t_vec)];
for ii = 1:numel(idx_set)
    % Grab Position
    this_idx = idx_set(ii);

    this_x = x_aoa(:,this_idx);
    this_v = v_aoa(:,this_idx);
    this_bng = atan2(this_v(2), this_v(1));

    % Make a triangular patch
    marker_radius = 250;
    num_pts = 3;
    vertex_theta = pi+(0:2*pi/num_pts:2*pi*(num_pts-1)/num_pts)-this_bng;
    marker_x = this_x(1) + marker_radius*cos(vertex_theta);
    marker_y = this_x(2) + marker_radius*sin(vertex_theta);

    % Draw Icon
    hdl = patch(marker_x, marker_y,'^','EdgeColor','k','FaceColor',hdl_traj.Color);
    utils.excludeFromLegend(hdl);

    % Draw LOB with uncertainty
    psi = triang.measurement(this_x, x_tgt);
    psi_high = psi + theta_unc * pi/180;
    psi_low =  psi - theta_unc * pi/180;

    xy_lob = triang.drawLob(this_x, psi, x_tgt, 5);
    xy_lob_high = triang.drawLob(this_x, psi_high, x_tgt, 5);
    xy_lob_low = triang.drawLob(this_x, psi_low, x_tgt, 5);
    lobFill = cat(2,xy_lob_high,fliplr(xy_lob_low),xy_lob_high(:,1));

    hdl_fill = fill(lobFill(1,:),lobFill(2,:),.2,'FaceAlpha',.2,'EdgeColor','k','LineStyle','--','DisplayName','Uncertainty Interval');
    hdl_lob = plot(xy_lob(1,:), xy_lob(2,:), '-.','DisplayName','LOB','Color',hdl_traj.Color);
    if ii > 1
        utils.excludeFromLegend(hdl_fill);
        utils.excludeFromLegend(hdl_lob);
    end
end
plot(x_tgt(1), x_tgt(2),'o','DisplayName','Target');
legend('Location','SouthEast');
ylim([-1e3 11e3]);
xlim([-7.5e3 13.5e3]);
utils.setPlotStyle(gca,{'widescreen'})

% Model CEP over time

% Since we can't do geolocation with the first measurement, let's
% initialize the track manually
x_prev = [0; 1e3];
P_prev = diag([1e3, 10e3]).^2;
C_df = (theta_unc*pi/180)^2;

cep_vec = zeros(size(t_vec));
for idx = 1:numel(t_vec)
    this_x_aoa = x_aoa(:,idx);
    
    z_fun = @(x) triang.measurement(this_x_aoa,x);
    H_fun = @(x) triang.jacobian(this_x_aoa,x)';

    this_psi = triang.measurement(this_x_aoa,x_tgt) + sqrt(C_df)*randn(1);

    [this_x, this_P] = tracker.ekfUpdate(x_prev, P_prev, this_psi, C_df, z_fun, H_fun);
    cep_vec(idx) = utils.computeCEP50(this_P);

    x_prev = this_x;
    P_prev = this_P;
end

fig8b=figure;
plot(t_vec, cep_vec/1e3);
xlabel('Time [s]');
ylabel('$CEP_{50}$ [m]');
utils.setPlotStyle(gca,{'widescreen'});

utils.exportPlot(fig8a,[prefix '8a']);
utils.exportPlot(fig8b,[prefix '8b']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;