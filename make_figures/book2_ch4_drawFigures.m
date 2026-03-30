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

%% Figures 10 and 11, Example 4.1
if force_recalc
    fprintf('Executing Example 4.1...\n');
    figs = book2_ex4_1();
    
    fig10=figs(1);
    figure(fig10);
    utils.setPlotStyle(gca,{'widescreen', 'tight'});
    utils.exportPlot(fig10,[prefix '10']);
    
    fig11=figs(2);
    figure(fig11);
    utils.setPlotStyle(gca,{'tight'});
    utils.exportPlot(fig11,[prefix '11'],true);
else
    fprintf('Skipping Example 4.1 (set force_recalc=true to run this example)...\n');
end

%% Figure 12, Example 4.2 
fprintf('Executing Example 4.2...\n');
fig12 = book2_ex4_2();

utils.setPlotStyle(gca,{'tight','equal'});
utils.exportPlot(fig12,[prefix '12'],true);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;