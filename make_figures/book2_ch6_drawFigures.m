% Draw Figures - Chapter 6
%
% This script generates all of the figures that appear in
% Chapter 6 of the textbook.
%
% Nicholas O'Donoughue
% 13 October 2021

% Clear Figures
close all;

% Set up directory and filename for figures
dirNm = fullfile(pwd,'figures','practical_geo');
if ~exist(dirNm,'dir')
    mkdir(dirNm);
end
prefix = fullfile(dirNm,'fig6_');

% Initialize Plot Preference
doFullColor=true;
utils.initPlotSettings;
colors=get(groot,'DefaultAxesColorOrder');

% Flag to force re-execution of long scripts
if ~exist('force_recalc','var')
    force_recalc = false;
end

% Reset the random number generator, to ensure reproducability
rng('default') ; 

addpath('examples');

%% Figure 6.1, Example 6.1
fig = book2_ex6_1;

utils.exportPlot(fig, [prefix '1']);

%% Figure 6.2, Example 6.2
fig = book2_ex6_2;

utils.exportPlot(fig, [prefix '2']);

%% Figure 6.2a, Example 6.3
fig = book2_ex6_3;

utils.exportPlot(fig, [prefix '2a']);

%% Figures 6.3a and 6.3b, Impact of Sensor Position Errors

x_aoa = [-1, 1;0, 0];
x_target = [0; 5];
x_err = [1;0];

x_aoa_un = x_aoa + x_err;
x_aoa_nonun = x_aoa + x_err*[-1 1];

% Solve True LOBs and AOAs
lob = x_target - x_aoa;
psi = triang.measurement(x_aoa, x_target);
r = utils.rng(x_aoa, x_target);

% Solve Erroneous AOAs
x_tgt_un = utils.find_intersect(x_aoa_un(:,1),psi(1),x_aoa_un(:,2),psi(2));
x_tgt_nonun = utils.find_intersect(x_aoa_nonun(:,1),psi(1),x_aoa_nonun(:,2),psi(2));

% num_dim x num_tx x 2
lob_zero = [cos(psi), sin(psi)]'.*reshape([zeros(2,1),3*r],1,2,2);

lob_true = x_aoa + lob_zero;
lob_un = x_aoa_un + lob_zero;
lob_nonun = x_aoa_nonun + lob_zero;

% Draw Uniform Offset case
fig3a=figure;
hdls = plot(squeeze(lob_true(1,:,:))',...
            squeeze(lob_true(2,:,:))','-','Color',colors(1,:),...
            'DisplayName','True LOB');
utils.excludeFromLegend(hdls(2:end));
hold on;
hdls = plot(squeeze(lob_un(1,:,:))',...
            squeeze(lob_un(2,:,:))','--','Color',colors(2,:),...
            'DisplayName','Perceived LOB');
utils.excludeFromLegend(hdls(2:end));
plot(x_aoa(1,:),x_aoa(2,:),'^','Color',colors(1,:),'DisplayName','True Sensor Position');
plot(x_aoa_un(1,:),x_aoa_un(2,:),'v','Color',colors(2,:),'DisplayName','Est. Sensor Position');
plot(x_target(1),x_target(2),'o','Color',colors(1,:),'DisplayName','Target');
plot(x_tgt_un(1),x_tgt_un(2),'o','Color',colors(2,:),'DisplayName','Est. Target');
xlim([-6 6]);
ylim([0 11]);
legend('Location','NorthWest');

utils.setPlotStyle(gca,{'tight','clean','equal'});
utils.exportPlot(fig3a,[prefix '3a']);

% Draw Non-Uniform Offset case
fig3b=figure;
plot(squeeze(lob_true(1,:,:))',...
     squeeze(lob_true(2,:,:))','-','Color',colors(1,:),'DisplayName','True LOB');
hold on;
plot(squeeze(lob_nonun(1,:,:))',...
     squeeze(lob_nonun(2,:,:))','--','Color',colors(2,:),'DisplayName','Perceived LOB');
plot(x_aoa(1,:),x_aoa(2,:),'^','Color',colors(1,:),'DisplayName','True Sensor Position');
plot(x_aoa_nonun(1,:),x_aoa_nonun(2,:),'v','Color',colors(2,:),'DisplayName','Est. Sensor Position');
plot(x_target(1),x_target(2),'o','Color',colors(1,:),'DisplayName','Target');
plot(x_tgt_nonun(1),x_tgt_nonun(2),'o','Color',colors(2,:),'DisplayName','Est. Target');
% legend('Location','NorthWest');
xlim([-6 6]);
ylim([0 11]);

utils.setPlotStyle(gca,{'tight','clean','equal'});
utils.exportPlot(fig3b,[prefix '3b']);

%% Figures 6.7a, 6.7b, and 6.8, Example 6.4
if force_recalc
    figs = book2_ex6_4;

    utils.exportPlot(figs(1), [prefix '7a']);
    utils.exportPlot(figs(2), [prefix '7b']);
    utils.exportPlot(figs(3), [prefix '8a']);
    utils.exportPlot(figs(4), [prefix '8b']);
end

%% Figure 6.10, Example 6.5
figs = book2_ex6_5;

utils.exportPlot(figs(1), [prefix '10']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;