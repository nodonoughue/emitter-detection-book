% Draw Figures - Chapter 5
%
% This script generates all of the figures that appear in
% Chapter 5 of the textbook.
%
% Nicholas O'Donoughue
% 1 July 2019

% Clear Figures
close all;

% Set up directory and filename for figures
dirNm = fullfile(pwd,'figures');
dirNmComponents = fullfile(pwd,'Graphics','Components');
if ~exist(dirNm,'dir')
    mkdir(dirNm);
end
if ~exist(dirNmComponents,'dir')
    mkdir(dirNmComponents);
end

prefix = fullfile(dirNm,'fig5_');

% Initialize Plot Preference
utils.initPlotSettings;

% Initializes colorSet - Mx3 RGB vector for successive plot lines
colors = get(0,'DefaultAxesColorOrder');

% Reset the random number generator, to ensure reproducability
rng('default') ; 

% Add path for folder with examples from the textbook
addpath('examples');

%% Figure 4 - Superhet Performance
fig4 = ex5_1;
utils.exportPlot(fig4,[prefix '4']);

%% Figure 6 - FMCW Radar
fig6 = ex5_2;
utils.exportPlot(fig6,[prefix '6']);

%% Figure 7 - Pulsed Radar
fig7 = ex5_3;
utils.exportPlot(fig7,[prefix '7']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;