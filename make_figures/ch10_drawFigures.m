% Draw Figures - Chapter 10
%
% This script generates all of the figures that appear in
% Chapter 10 of the textbook.
%
% Nicholas O'Donoughue
% 1 July 2019

% Clear Figures
close all;

% Flag to force re-execution of long scripts
if ~exist('force_recalc','var')
    force_recalc = false;
end

% Set up directory and filename for figures
dirNm = fullfile(pwd,'figures');
if ~exist(dirNm,'dir')
    mkdir(dirNm);
end
prefix = fullfile(dirNm,'fig10_');

% Initialize Plot Preference
utils.initPlotSettings;

% Reset the random number generator, to ensure reproducability
rng('default') ; 

% Add folder for examples from the textbook
addpath('examples');

%% Figure 1, LOB Intersection
% Compute an isochrone between sensors 1 and 2, then draw an AOA slice from
% sensor 3

% Initialize Detector/Source Locations
x0 = [0; 0];
x1 = [1;1];
x2 = [-1; 0];

x_source = [.1; .9];

% Compute Ranges
r0 = utils.rng(x0,x_source);
r1 = utils.rng(x1,x_source);
r2 = utils.rng(x2,x_source);

% Error Values
epsang = 5*pi/180;

% Find AOA 
lob0 = x_source - x0;
aoa0 = atan2(lob0(2),lob0(1));
xaoa0 = x0 + [0 cos(aoa0);0 sin(aoa0)]*5*r1;
xaoap0 = x0 + [0 cos(aoa0+epsang);0 sin(aoa0+epsang)]*5*r1;
xaoam0 = x0 + [0 cos(aoa0-epsang);0 sin(aoa0-epsang)]*5*r1;
lobFill0 = cat(2,xaoap0,fliplr(xaoam0),xaoap0(:,1));

lob1 = x_source - x1;
aoa1 = atan2(lob1(2),lob1(1));
xaoa1 = x1 + [0 cos(aoa1);0 sin(aoa1)]*5*r1;
xaoap1 = x1 + [0 cos(aoa1+epsang);0 sin(aoa1+epsang)]*5*r1;
xaoam1 = x1 + [0 cos(aoa1-epsang);0 sin(aoa1-epsang)]*5*r1;
lobFill1 = cat(2,xaoap1,fliplr(xaoam1),xaoap1(:,1));

lob2 = x_source - x2 + [.1;0];

aoa2 = atan2(lob2(2),lob2(1));
xaoa2 = x2 + [0 cos(aoa2);0 sin(aoa2)]*5*r2;
xaoap2 = x2 + [0 cos(aoa2+epsang);0 sin(aoa2+epsang)]*5*r2;
xaoam2 = x2 + [0 cos(aoa2-epsang);0 sin(aoa2-epsang)]*5*r2;
lobFill2 = cat(2,xaoap2,fliplr(xaoam2),xaoap2(:,1));

% Draw Figure
fig1 = figure();hold on;

% LOBs
plot(xaoa0(1,:),xaoa0(2,:),'k-','DisplayName','AOA Solution');
h=plot(xaoa1(1,:),xaoa1(2,:),'k-');
h1=plot(xaoa2(1,:),xaoa2(2,:),'k-');
utils.excludeFromLegend([h,h1]);

% Uncertainty Intervals
plot(xaoap0(1,:),xaoap0(2,:)','k--','DisplayName','Uncertainty Interval');
h=plot(xaoam0(1,:),xaoam0(2,:),'k--');
utils.excludeFromLegend(h);
h = fill(lobFill0(1,:),lobFill0(2,:),'k','FaceAlpha',.1,'EdgeColor','none');
utils.excludeFromLegend(h);

h=plot(xaoap1(1,:),xaoap1(2,:)','k--',xaoam1(1,:),xaoam1(2,:),'k--');
utils.excludeFromLegend(h);
h = fill(lobFill1(1,:),lobFill1(2,:),'k','FaceAlpha',.1,'EdgeColor','none');
utils.excludeFromLegend(h);

h=plot(xaoap2(1,:),xaoap2(2,:)','k--',xaoam2(1,:),xaoam2(2,:),'k--');
utils.excludeFromLegend(h);
h = fill(lobFill2(1,:),lobFill2(2,:),'k','FaceAlpha',.1,'EdgeColor','none');
utils.excludeFromLegend(h);


% Position Markers
plot([x0(1),x1(1),x2(1)],[x0(2),x1(2),x2(2)],'ko','DisplayName','Sensors');
%plot(xs(1),xs(2),'k^','MarkerSize',8,'DisplayName','Transmitter');

% Position Labels
text(x0(1)+.05,x0(2)-.1,'$S_0$');
text(x1(1)+.05,x1(2)-.1,'$S_1$');
text(x2(1)+.05,x2(2)-.1,'$S_2$');

% Adjust Axes
legend('Location','SouthEast');
utils.setPlotStyle(gca,{'clean','equal','widescreen','tight'});
ylim([-.5 1.5]);
xlim([-1 2]);

utils.exportPlot(fig1,[prefix '1']);

%% Figure 2 - This time with angles

% Initialize Detector/Source Locations
x0 = [0; 0];
x1 = [1;1];
x2 = [-1; 0];
x_source = [.1; .9];

% Compute Ranges
r0 = utils.rng(x0,x_source);
r1 = utils.rng(x1,x_source);
r2 = utils.rng(x2,x_source);

% Find AOA 
lob0 = x_source - x0;
aoa0 = atan2(lob0(2),lob0(1));
xaoa0 = x0 + [0 cos(aoa0);0 sin(aoa0)]*5*r1;

lob1 = x_source - x1;
aoa1 = atan2(lob1(2),lob1(1))+2*pi;
xaoa1 = x1 + [0 cos(aoa1);0 sin(aoa1)]*5*r1;

lob2 = x_source - x2;
aoa2 = atan2(lob2(2),lob2(1));
xaoa2 = x2 + [0 cos(aoa2);0 sin(aoa2)]*5*r2;

% Draw Figure
fig2 = figure();hold on;

% LOBs
h=plot(xaoa0(1,:),xaoa0(2,:),'k-','DisplayName','AOA Solution');
utils.excludeFromLegend(h);
h=plot(xaoa1(1,:),xaoa1(2,:),'k-');
utils.excludeFromLegend(h);
h=plot(xaoa2(1,:),xaoa2(2,:),'k-');
utils.excludeFromLegend(h);

% Angle Markers
angle_rad = .1;
h0=plot(x0(1)+[0 3*angle_rad],x0(2)*[1 1],'k-','LineWidth',.5);
h1=plot(x0(1)+angle_rad*cos(linspace(0,aoa0,100)),...
        x0(2)+angle_rad*sin(linspace(0,aoa0,100)),'k-','LineWidth',.5);
text(x0(1)+angle_rad,x0(2)+angle_rad,'$\psi_0$');
h2=plot(x1(1)+[0 3*angle_rad],x1(2)*[1 1],'k-','LineWidth',.5);
h3=plot(x1(1)+angle_rad*cos(linspace(0,aoa1,100)),...
        x1(2)+angle_rad*sin(linspace(0,aoa1,100)),'k-','LineWidth',.5);
text(x1(1)+angle_rad,x1(2)+angle_rad,'$\psi_1$');
h4=plot(x2(1)+[0 3*angle_rad],x2(2)*[1 1],'k-','LineWidth',.5);
h5=plot(x2(1)+angle_rad*cos(linspace(0,aoa2,100)),...
        x2(2)+angle_rad*sin(linspace(0,aoa2,100)),'k-','LineWidth',.5);
text(x2(1)+2*angle_rad,x2(2)+angle_rad,'$\psi_2$');
utils.excludeFromLegend([h0,h1,h2,h3,h4,h5]);

% Position Markers
plot([x0(1),x1(1),x2(1)],[x0(2),x1(2),x2(2)],'ko','DisplayName','Sensors');
plot(x_source(1),x_source(2),'k^','MarkerSize',8,'DisplayName','Transmitter');

% Position Labels
text(x0(1)+.05,x0(2)-.1,'$(x_0,y_0)$');
text(x1(1)+.05,x1(2)-.1,'$(x_1,y_1)$');
text(x2(1)+.05,x2(2)-.1,'$(x_2,y_2)$');
text(x_source(1)+.05,x_source(2)-.1,'$(x,y)$');

% Adjust Axes
legend('Location','SouthEast');
utils.setPlotStyle(gca,{'clean','equal','widescreen','tight'});
ylim([-.5 1.5]);
xlim([-1 2]);

utils.exportPlot(fig2,[prefix '2']);

%% Figure 3 - Geometric Solutions

% Define triangle
x_sensor = [0 1 .5];
y = [0 1 1.5];

% Initialize all 3 plots
fig3=figure;
for idx_plot=1:2
    subplot(1,2,idx_plot);
    plot([x_sensor(1:end),x_sensor(1)],[y(1:end),y(1)]);
    axis equal
    hold on;
end

% Midpoint - connect each midoint to opposite vertex
mx = .5 * (x_sensor + circshift(x_sensor,-1,2));
my = .5 * (y + circshift(y,-1,2));

% Vector to opposite vertex
v = [circshift(x_sensor,1,2);circshift(y,1,2)] - [mx;my];

subplot(1,2,1);
for i=1:numel(mx)
    plot(mx(i) + [0 v(1,i)],my(i)+[0 v(2,i)],'k-','LineWidth',.5);
end

% Angle Bisector method
th_fwd = atan2(circshift(y,-1,2)-y,circshift(x_sensor,-1,2)-x_sensor);
th_back = atan2(circshift(y,1,2)-y,circshift(x_sensor,1,2)-x_sensor);
th_bisector = .5*(th_fwd+th_back);
flipped = th_bisector < th_fwd;
th_bisector(flipped) = th_bisector(flipped)+pi;

% Find Intersection on Opposite Edge
R_fwd = sqrt(abs(x_sensor - circshift(x_sensor,-1,2)).^2+...
             abs(y - circshift(y,-1,2)).^2);
R_opp = circshift(R_fwd,-1,2);
R_rev = circshift(R_fwd,1,2);

R_bisect = sqrt((R_fwd.*R_rev).*((R_fwd+R_rev).^2-R_opp.^2)./(R_fwd+R_rev).^2);

% Plot
subplot(1,2,2);
for i=1:numel(R_bisect)
    plot(x_sensor(i)+ R_bisect(i)*[0 cos(th_bisector(i))],...
         y(i) + R_bisect(i)*[0 sin(th_bisector(i))],'k-','LineWidth',.5);
end

for i=1:2
    subplot(1,2,i);
    utils.setPlotStyle(gca,{'equal','clean'});
end

utils.exportPlot(fig3,[prefix '3']);

%% Figure 4 - Example Results
% Test LOB intersection

x_sensor= [0 2;
    1 0;
    2 0]';

psi = [65 75 95]*pi/180;
r = 0:20;

% Plot the points and lobs
fig4=figure;
hh=plot(x_sensor(1,:),x_sensor(2,:),'ko','DisplayName','Sensors');
utils.excludeFromLegend(hh);
text(x_sensor(1,1)-.4,x_sensor(2,1),'$S_0$');
text(x_sensor(1,2)-.4,x_sensor(2,2),'$S_1$');
text(x_sensor(1,3)-.4,x_sensor(2,3),'$S_2$');
hold on;
h=plot((x_sensor(1,:)' + cos(psi(:))*r)',...
       (x_sensor(2,:)' + sin(psi(:))*r)','k-','LineWidth',.5,'DisplayName','AOA Solutions');
utils.excludeFromLegend(h);

% Look for intersections
x_centroid = triang.centroid(x_sensor,psi);
x_angleBisector = triang.angle_bisector(x_sensor,psi);

plot(x_centroid(1),x_centroid(2),'kx','DisplayName','Centroid');
plot(x_angleBisector(1),x_angleBisector(2),'k+','DisplayName','Incenter');
legend('Location','NorthWest');
xlim([-1 5]);
ylim([0 10]);
utils.setPlotStyle(gca,{'clean','widescreen','tight'});
utils.exportPlot(fig4,[prefix '4']);

%% Figure 5 - Example Solution
if force_recalc
    
[fig5a,fig5b]=ex10_1();

utils.exportPlot(fig5a,[prefix,'5a']);
utils.exportPlot(fig5b,[prefix,'5b']);

end

%% Figure 6 - 2 Sensor Configuration
fig6 = ex10_2;
utils.exportPlot(fig6,[prefix '6']);

%% Figure 7 - 3 Sensor Configuration
fig7 = ex10_3;
utils.exportPlot(fig7,[prefix '7']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;