% Draw Figures - Chapter 6
%
% This script generates all of the figures that appear in
% Chapter 6 of the textbook.
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

prefix = fullfile(dirNm,'fig6_');

% Initialize Plot Preference
utils.initPlotSettings;

% Initializes colorSet - Mx3 RGB vector for successive plot lines
colors = get(0,'DefaultAxesColorOrder');

% Reset the random number generator, to ensure reproducability
rng('default') ; 

%% Figure 1, Bayesian Example
d_lam = 4;
N = 10;
psi_0 = 135*pi/180;
%psi_0 = 95*pi/180;

% Narrow Beam (Marginal Distribution)
f_el = @(psi) abs(cos(psi-pi/2)).^1.2;
f_af = @(psi) abs(sin(pi*d_lam*N*(cos(psi)-cos(psi_0)))./sin(pi*d_lam*(cos(psi)-cos(psi_0))))/N;

% Wide Beam (Prior Distribution)
d_lam = .5;
f_af_2 = @(psi) abs(sin(pi*d_lam*N*(cos(psi)-cos(psi_0)))./sin(pi*d_lam*(cos(psi)-cos(psi_0))))/N;

fig1 =figure;
psi_vec = 0:pi/1001:pi;
plot(psi_vec,f_el(psi_vec).*f_af(psi_vec),'DisplayName','Narrow Beam (marginal)');
hold on;set(gca,'ColorOrderIndex',4);
plot(psi_vec,f_af_2(psi_vec),'DisplayName','Wide Beam (prior)');
grid on;
xlabel('$\psi$');
legend('Location','NorthWest');
utils.setPlotStyle(gca,{'widescreen'});
utils.exportPlot(fig1,[prefix '1']);

%% Figure 2, Convex Optimization Example

x0 = [1,.5];

xx = -5:.01:5;
yy = -3:.01:3;

[XX,YY] = ndgrid(xx,yy);

f = 1.*(XX-x0(1)).^2 + 5.*(YY-x0(2)).^2;

% Initial Estimate
x_est = [-3, 2;
          2,-1.5;
          1, 2.2;
          1,  .6];

fig2=figure;
contour(XX,YY,f);colormap(gray);
hold on;
plot(x0(1),x0(2),'^');
plot(x_est(:,1),x_est(:,2),'--+');

text(-3,2.1,'Initial Estimate','FontSize',10);
text(-.4,.25,'True Minimum','FontSize',10);

utils.setPlotStyle(gca,{'clean','widescreen','tight'});
utils.exportPlot(fig2,[prefix '2']);

%% Figure 3, Tracker Example

y = [1   1.1  1.3 1.4  1.35  1.3   .7 .75];
x = [1   1.05 1.2 1.35 1.45 1.35 1.2 .8];

s2 = [.8  .5   .4  .3   .3   .2   .2  .6];
N = numel(y);

fig3=figure;
hold on;
for i=2:numel(s2)
    h=fill(i + [.2 .2 -.2 -.2 .2],...
          x(i-1) + s2(i)*[-1 1 1 -1 -1],...
          [.8,.8,.8],'EdgeColor','none');
    utils.excludeFromLegend(h);
end
plot(1:N,y,'kx','DisplayName','Measurement');
plot(1:N,x,'ko-.','DisplayName','Estimate');
legend('Location','NorthWest');
xlabel('Time');
ylabel('Parameter ($\theta$)');

utils.setPlotStyle(gca,{'notick','widescreen'});
utils.exportPlot(fig3,[prefix '3']);

%% Figure 4, Angle Error Variance

% Sensor Coordinates
x0 = [0;0];
xs = [2,1.4];
% Bearing
psi_deg = 40;
std_dev = 8;
len=5;
aoa_lower = psi_deg - std_dev;
aoa_upper = psi_deg + std_dev;

fig4=figure;
plot(x0(1),x0(2),'o');
hold on;
fill(x0(1) + [0 len*cosd(aoa_lower) len*cosd(psi_deg) len*cosd(aoa_upper) 0],...
     x0(2) + [0 len*sind(aoa_lower) len*sind(psi_deg) len*sind(aoa_upper) 0],[.8,.8,.8]);
 
plot(x0(1)+[0 len*cosd(psi_deg)],x0(2)+[0 len*sind(psi_deg)],'k-.');
plot(x0(1) +[0 len*cosd(aoa_lower)],x0(2)+[0 len*sind(aoa_lower)],'k-');
plot(x0(1) +[0 len*cosd(aoa_upper)],x0(2)+[0 len*sind(aoa_upper)],'k-');

text(-.15,.1,'Sensor','FontSize',10);
text(.95,.85,'Estimated Bearing','FontSize',10,'Rotation',35);
%plot([1 1.2],[1.4,1],'k-','LineWidth',.5);
text(1.6,1.6,'Confidence Interval','FontSize',10,'BackgroundColor',[.8,.8,.8]);
%plot([1.45 1.6],[.65 1.1],'k-','LineWidth',.5);

plot(len/5*cosd(linspace(psi_deg,aoa_upper,10)),len/5*sind(linspace(psi_deg,aoa_upper,10)),'k-','LineWidth',.5);
text(len/5*cosd(aoa_upper)-.2,len/5*sind(aoa_upper)+.1,'RMSE','FontSize',10);

%plot(xs(1),xs(2),'x');
%text(xs(1)-.1,xs(2)+.1,'Source','FontSize',10);

xlim([-.5 2.5]);
ylim([-.2 1.8]);
utils.setPlotStyle(gca,{'clean','widescreen','tight'});
utils.exportPlot(fig4,[prefix '4']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;