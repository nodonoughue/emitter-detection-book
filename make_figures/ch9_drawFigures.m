% Draw Figures - Chapter 9
%
% This script generates all of the figures that appear in
% Chapter 9 of the textbook.
%
% Nicholas O'Donoughue
% 1 July 2019

% Clear Figures
close all;

% Set up directory and filename for figures
dirNm = fullfile(pwd,'figures');
if ~exist(dirNm,'dir')
    mkdir(dirNm);
end
prefix = fullfile(dirNm,'fig9_');

% Initialize Plot Preference
utils.initPlotSettings;

% Add path to folder for examples from the textbook
addpath('examples');

%% Figure 1, Plot of Error Ellipse

% Define Positions
x0 = [0;1];
x2 = [.2,;.2];

% Define Covariance Matrix
sx2 = 5;
sy2 = 3;
rho = .8;
sxy = rho*sqrt(sx2*sy2); % cross-covariance
C = [sx2 sxy; sxy sy2];

% Compute Elippses
xellipse1 = utils.drawErrorEllipse(x2,C,361,1); % 1 sigma
xellipse2 = utils.drawErrorEllipse(x2,C,361,95);% 95% confidence interval

% Draw figure
fig1 = figure();hold on;

% Draw True/Estimate Positions
plot(x0(1),x0(2),'k^','DisplayName','True');
plot(x2(1),x2(2),'k+','DisplayName','Estimated');

% Draw error ellipses
plot(xellipse1(1,:),xellipse1(2,:),'k-','DisplayName','1\sigma Ellipse');
plot(xellipse2(1,:),xellipse2(2,:),'k--','DisplayName','95% Ellipse');

% Adjust Figure Display
xlim([-1 1]);
%ylim([0 2]);
legend('Location','NorthWest');
utils.setPlotStyle(gca,{'clean','equal','widescreen','tight'});

% Save Figure
utils.exportPlot(fig1,[prefix '1']);

%% Figure 2 Plot of Error Ellipse Example
fig2=ex9_1;
utils.exportPlot(fig2,[prefix '2']);

%% Figure 3, Plot of CEP50

% Initialize Emitter Location and Estimate
x0 = [0;0];
x2 = [.5;-.2];

% Initialize Covariance Matrix
sx2 = 5;
sy2 = 3;
rho = .8;
sxy = rho*sqrt(sx2*sy2);
C = [sx2 sxy; sxy sy2];

% Compute Error Ellipses
xellipse = utils.drawErrorEllipse(x2,C,361,50);
xcep = utils.drawCEP50(x2,C,361);

% Draw Figure
fig3 = figure();hold on;

% Draw Ellipses
plot(xcep(1,:),xcep(2,:),'k-','DisplayName','CEP_{50}')
plot(xellipse(1,:),xellipse(2,:),'k--','DisplayName','50% Error Ellipse')
plot(x0(1),x0(2),'k^','DisplayName','True')
plot(x2(1),x2(2),'k+','DisplayName','Estimated')
xlim(1.1*max(xellipse(1,:))*[-1 1]);

% Adjust Display
legend('Location','NorthWest');
utils.setPlotStyle(gca,{'clean','equal','widescreen','tight'});

utils.exportPlot(fig3,[prefix '3']);

%% Figure 4 Plot of CEP50 and Error Ellipse Example
fig4=ex9_2;
utils.exportPlot(fig4,[prefix '4']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;