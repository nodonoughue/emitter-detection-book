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

%% Figure 1, Logarithmic Barriers

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
utils.exportPlot(fig1,[prefix '1']);


%% Illustration of Likelihood and Priors

% Set up axes
y = linspace(0,10,101);
x = linspace(-10,10,1001);

[xx,yy] = meshgrid(x,y);

% Data Likelihood
y_ctr = 5;
x_ctr = 0;
mean_vec = x_ctr - abs(y-y_ctr).^2/4;
mean = x_ctr - abs(yy-y_ctr).^2/4;
std_dev = 3;

ell = normpdf(xx,mean,std_dev);

fig2a=figure;
imagesc(x,y,ell)
hold on;
plot(mean_vec, y, 'k--');
grid on;

xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('Data Likelihood');

utils.setPlotStyle(gca,{'widescreen'});
utils.exportPlot(fig2a,[prefix '2a']);

% A priori pdf

prior = reshape(mvnpdf([xx(:), yy(:)], [x_ctr, y_ctr], [20, 2]),size(xx));

fig2b=figure;
imagesc(x,y,prior);
grid on;

xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('Prior Distribution on Target Location');

utils.setPlotStyle(gca,{'widescreen'});
utils.exportPlot(fig2b,[prefix '2b']);

% Posterier
fig2c=figure;
imagesc(x,y,ell.*prior);
grid on;

xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('Posterior Distribution on Target Location');

utils.setPlotStyle(gca,{'widescreen'});
utils.exportPlot(fig2c,[prefix '2c']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;