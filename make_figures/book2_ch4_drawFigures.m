% Draw Figures - Chapter 4
%
% This script generates all of the figures that appear in
% Chapter 4 of the textbook.
%
% Nicholas O'Donoughue
% 28 August 2021

% Clear Figures
close all;

% Set up directory and filename for figures
dirNm = fullfile(pwd,'figures','practical_geo');
if ~exist(dirNm,'dir')
    mkdir(dirNm);
end
prefix = fullfile(dirNm,'fig4_');

% Initialize Plot Preference
doFullColor=true;
utils.initPlotSettings;

% Reset the random number generator, to ensure reproducability
rng('default') ; 

addpath('examples');

%% Figures 9 and 10, Example 4.1
if force_recalc
    figs = book2_ex4_1();
    
    fig9=figs(1);
    figure(fig9);
    utils.setPlotStyle(gca,{'widescreen', 'tight'});
    utils.exportPlot(fig9,[prefix '9']);
    
    fig10=figs(2);
    figure(fig10);
    utils.setPlotStyle(gca,{'tight'});
    utils.exportPlot(fig10,[prefix '10']);
end

%% Figures 11 and 12, Example 4.2 
fig11 = book2_ex4_2();

utils.setPlotStyle(gca,{'tight','equal'});
utils.exportPlot(fig11,[prefix '11']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;