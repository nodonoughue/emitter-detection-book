% Draw Figures - Chapter 8
%
% This script generates all of the figures that appear in
% Chapter 8 of the textbook.
%
% Nicholas O'Donoughue
% 25 January 2022

% Clear Figures
close all;

% Set up directory and filename for figures
dirNm = fullfile(pwd,'figures','practical_geo');
if ~exist(dirNm,'dir')
    mkdir(dirNm);
end
prefix = fullfile(dirNm,'fig8_');

% Initialize Plot Preference
doFullColor=true;
utils.initPlotSettings;
colors=get(groot,'DefaultAxesColorOrder');

% Reset the random number generator, to ensure reproducability
rng('default'); 

addpath('examples');

%% Figures 8.3 and 8.4, Example 8.1
figs = book2_ex8_1;

utils.exportPlot(figs(1), [prefix '3']);
utils.exportPlot(figs(2), [prefix '4']);
utils.exportPlot(figs(3), [prefix '5']);

%% Figures 8.6a/b and 8.7, Example 8.2

figs = book2_ex8_2;

utils.exportPlot(figs(1), [prefix '6a']);
utils.exportPlot(figs(2), [prefix '6b']);
utils.exportPlot(figs(3), [prefix '7']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;