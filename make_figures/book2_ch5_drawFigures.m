    % Draw Figures - Chapter 5
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
prefix = fullfile(dirNm,'fig5_');

% Initialize Plot Preference
doFullColor=true;
utils.initPlotSettings;

% Reset the random number generator, to ensure reproducability
rng('default') ; 

addpath('examples');

%% Figures 5 and 6, Example 5.1
fprintf('Executing Example 5.1..\n');
figs = book2_ex5_1;

utils.exportPlot(figs(1),[prefix '5']);
utils.exportPlot(figs(2),[prefix '6']);
utils.exportPlot(figs(3),[prefix '6b']);

%% Figure 8, Example 5.2
fprintf('Executing Example 5.2..\n');
figs = book2_ex5_2;

utils.exportPlot(figs(1), [prefix '8']);
utils.exportPlot(figs(2), [prefix '8b']);

%% Figure 10, Example 5.3
fprintf('Executing Example 5.3..\n');
fig = book2_ex5_3;

utils.exportPlot(fig, [prefix '10']);

%% Figure 13, Example 5.4
fprintf('Executing Example 5.4..\n');
fig = book2_ex5_4;

utils.exportPlot(fig, [prefix '13']);

%% Figure 14, Logarithmic Barriers

% Ideal barrier
u = linspace(-3,1,1001);
I = zeros(size(u));
I(u>0) = inf;
I(abs(u)<=.001) = 10;

% Log barriers
t = [.5, 1, 2];
l = -(1./t(:)).*log(u);
l(:,u>0) = inf;

fig1=figure;
plot(u,I,'k--');
hold on;
plot(u, real(l));
legend('Ideal','t=0.5','t=1.0','t=1.5');
ylim([-5 5]);
xlim([-3 1])
xlabel('u');
ylabel('Cost');

utils.setPlotStyle(gca,{'widescreen', 'tight'});
utils.exportPlot(fig1,[prefix '14']);


%% Figure 15, Illustration of Likelihood and Priors

% Set up axes
y = linspace(0,10,101);
x = linspace(-10,10,1001);

[xx,yy] = meshgrid(x,y);

% Data Likelihood
y_ctr = 5;
x_ctr = 0;
mean_vec = x_ctr - abs(y-y_ctr).^2/4;
mean_val = x_ctr - abs(yy-y_ctr).^2/4;
std_dev = 3;

ell = normpdf(xx,mean_val,std_dev);

fig13a=figure;
imagesc(x,y,ell)
hold on;
plot(mean_vec, y, 'k--');
grid on;

xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('Data Likelihood');

utils.setPlotStyle(gca,{'widescreen'});
utils.exportPlot(fig13a,[prefix '15a']);

% A priori pdf

prior = reshape(utils.mvnpdf([xx(:), yy(:)], [x_ctr, y_ctr], [20, 2]),size(xx));

fig13b=figure;
imagesc(x,y,prior);
grid on;

xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('Prior Distribution on Target Location');

utils.setPlotStyle(gca,{'widescreen'});
utils.exportPlot(fig13b,[prefix '15b']);

% Posterier
fig13c=figure;
imagesc(x,y,ell.*prior);
grid on;

xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('Posterior Distribution on Target Location');

utils.setPlotStyle(gca,{'widescreen'});
utils.exportPlot(fig13c,[prefix '15c']);

%% Figure 16, Example 5.5
fprintf('Executing Example 5.5...\n');
figs = book2_ex5_5;

utils.exportPlot(figs(1), [prefix '16a']);
utils.exportPlot(figs(2), [prefix '16b']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;