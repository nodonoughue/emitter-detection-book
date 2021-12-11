% Draw Figures - Chapter 2
%
% This script generates all of the figures that appear in
% Chapter 2 of the textbook.
%
% Nicholas O'Donoughue
% 16 April 2021

% Clear Figures
close all;

% Set up directory and filename for figures
dirNm = fullfile(pwd,'figures','practical_geo');
if ~exist(dirNm,'dir')
    mkdir(dirNm);
end
prefix = fullfile(dirNm,'fig2_');

% Initialize Plot Preference
doFullColor=true;
utils.initPlotSettings;

% Reset the random number generator, to ensure reproducability
rng('default') ; 

if ~exist('force_recalc','var')
    force_recalc = false;
end

%% Figure 1, AOA Geometry (recreation of figure 10.2)

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
fig1 = figure();hold on;

% LOBs
h=plot(xaoa0(1,:),xaoa0(2,:),'k-','DisplayName','Line of Bearing');
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

utils.exportPlot(fig1,[prefix '1']);

%% Figure 2, TDOA Geometry (recreation of figure 11.1b)
% Initialize Detector/Source Locations
x_sensor1 = [0; 0];
x_sensor2 = [.8; .2];
x_sensor3 = [1;1];
x_source = [.1; .9];

% Compute Ranges
r1 = utils.rng(x_sensor1,x_source);
r2 = utils.rng(x_sensor2,x_source);
r3 = utils.rng(x_sensor3,x_source);

% Find Isochrones
xiso1 = tdoa.drawIsochrone(x_sensor1,x_sensor2,r2-r1,1000,3);
xiso2 = tdoa.drawIsochrone(x_sensor2,x_sensor3,r3-r2,1000,3);

% Draw Figure
fig2 = figure();hold on;

% Isochrones
plot(xiso1(1,:),xiso1(2,:),'k:','DisplayName','Isochrone');
hiso2=plot(xiso2(1,:),xiso2(2,:),'k:');
utils.excludeFromLegend(hiso2);

% Isochrone Labels
text(mean([x_sensor1(1),x_sensor2(1)]),mean([x_sensor1(2),x_sensor2(2)])-.2,'$TDOA_{1,2}$');
text(mean([x_sensor2(1),x_sensor3(1)])+.3,mean([x_sensor2(2),x_sensor3(2)]),'$TDOA_{2,3}$');

% Position Markers
hiso2=plot([x_sensor1(1),x_sensor2(1),x_sensor3(1)],[x_sensor1(2),x_sensor2(2),x_sensor3(2)],'k-','LineWidth',1);
utils.excludeFromLegend(hiso2);
plot([x_sensor1(1),x_sensor2(1),x_sensor3(1)],[x_sensor1(2),x_sensor2(2),x_sensor3(2)],'ko','DisplayName','Sensors');
plot(x_source(1),x_source(2),'k^','MarkerSize',8,'DisplayName','Transmitter');

% Position Labels
text(x_sensor1(1)+.05,x_sensor1(2)-.1,'$S_1$');
text(x_sensor2(1)+.05,x_sensor2(2)-.1,'$S_2$');
text(x_sensor3(1)+.05,x_sensor3(2)-.1,'$S_3$');

% text(x_source(1)+.05,x_source(2)+.05,'$T_1$');

% Adjust Axes
xlim([-2 3]);
ylim([-1 2]);
legend('Location','SouthEast');

utils.setPlotStyle(gca,{'clean','equal','widescreen','tight'});
utils.exportPlot(fig2,[prefix '2']);

%% Figure 3, FDOA Geometry (recreation of figure 12.1)
x_source = [0;2]; % Transmitter/source
x_sensor = [-1,0;
       0,0;
       2,0]';
v_sensor = [1,1;
      1,1;
      1,1]';


% Draw Geometry
fig3=figure;
plot(x_source(1),x_source(2),'k^','DisplayName','Transmitter');
hold on;
plot(x_sensor(1,:),x_sensor(2,:),'ko','DisplayName','Sensors');
text(x_sensor(1,1)-.2,x_sensor(2,1)-.2,'$S_1$','FontSize',10);
text(x_sensor(1,2)-.2,x_sensor(2,2)-.2,'$S_2$','FontSize',10);
text(x_sensor(1,3)-.2,x_sensor(2,3)-.2,'$S_3$','FontSize',10);

% Draw velocity arrows
utils.drawArrow(x_sensor(1,1)+[0 v_sensor(1,1)]/4,x_sensor(2,1)+[0 v_sensor(2,1)]/4);
utils.drawArrow(x_sensor(1,2)+[0 v_sensor(1,2)]/4,x_sensor(2,2)+[0 v_sensor(2,2)]/4);
utils.drawArrow(x_sensor(1,3)+[0 v_sensor(1,3)]/4,x_sensor(2,3)+[0 v_sensor(2,3)]/4);

% Draw isodoppler line S12
vdiff12 = utils.dopDiff(x_source,[0 0]',x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),3e8);
x_isodop12 = fdoa.drawIsodop(x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),vdiff12,1000,5);
plot(x_isodop12(1,:),x_isodop12(2,:),'k-.','DisplayName','Line of Constant FDOA');
text(-1,2.7,'$S_{12}$ Solution','FontSize',10);

% Draw isodoppler line S23
vdiff23 = utils.dopDiff(x_source,[0 0]',x_sensor(:,2),v_sensor(:,2),x_sensor(:,3),v_sensor(:,3),3e8);
x_isodop23 = fdoa.drawIsodop(x_sensor(:,2),v_sensor(:,2),x_sensor(:,3),v_sensor(:,3),vdiff23,1000,5);
hh=plot(x_isodop23(1,:),x_isodop23(2,:),'k-.');
utils.excludeFromLegend(hh);
text(1.5,.85,'$S_{23}$ Solution','FontSize',10);
xlim([-2 3]);
ylim([-1 4]);
legend;

utils.setPlotStyle(gca,{'clean','widescreen','tight'});
utils.exportPlot(fig3,[prefix '3']);


%% Figure 4, Hybrid Geometry (recreation of figure 13.1)
x_source = [2;4]; % Transmitter/source

x_sensor = [0,0; 3,0]'; % Receivers
v_sensor = [1,0;1,0]';

fig4=figure;
plot(x_source(1),x_source(2),'k^','DisplayName','Source');
hold on;

% Position Markers
plot(x_sensor(1,:),x_sensor(2,:),'ko','DisplayName','Sensors');

% Draw velocity arrows
utils.drawArrow(x_sensor(1,1)+[0 v_sensor(1,1)]/4,x_sensor(2,1)+[0 v_sensor(2,1)]/4);
utils.drawArrow(x_sensor(1,2)+[0 v_sensor(1,2)]/4,x_sensor(2,2)+[0 v_sensor(2,2)]/4);

% -- Direction of Arrival

% Compute Ranges
r = utils.rng(x_sensor,x_source);

% Error Values
epsang = 5*pi/180;

% Find AOA 
lob = x_source - x_sensor;
aoa = atan2(lob(2,:),lob(1,:));

xaoa1 = x_sensor(:,1) + [0 cos(aoa(1));0 sin(aoa(1))]*5*r(1);
xaoap1 = x_sensor(:,1) + [0 cos(aoa(1)+epsang);0 sin(aoa(1)+epsang)]*5*r(1);
xaoam1 = x_sensor(:,1) + [0 cos(aoa(1)-epsang);0 sin(aoa(1)-epsang)]*5*r(1);
lobFill1 = cat(2,xaoap1,fliplr(xaoam1),xaoap1(:,1));

xaoa2 = x_sensor(:,2) + [0 cos(aoa(2));0 sin(aoa(2))]*5*r(2);
xaoap2 = x_sensor(:,2) + [0 cos(aoa(2)+epsang);0 sin(aoa(2)+epsang)]*5*r(2);
xaoam2 = x_sensor(:,2) + [0 cos(aoa(2)-epsang);0 sin(aoa(2)-epsang)]*5*r(2);
lobFill2 = cat(2,xaoap2,fliplr(xaoam2),xaoap2(:,1));

% LOBs
plot(xaoa1(1,:),xaoa1(2,:),'k-','DisplayName','Line of Bearing');
h=plot(xaoa2(1,:),xaoa2(2,:),'k-');
utils.excludeFromLegend(h);

% -- Time Difference of Arrival
% Initialize Detector/Source Locations

% Isochrones
dr = utils.rngDiff(x_source,x_sensor(:,1),x_sensor(:,2));
xiso = tdoa.drawIsochrone(x_sensor(:,1),x_sensor(:,2),dr,1000,30);
plot(xiso(1,:),xiso(2,:),':','DisplayName','Line of Constant TDOA');
% text(.5,-.5,'TDOA Solution');

% -- Frequency Difference of Arrival

% Draw isodoppler line
vdiff = utils.dopDiff(x_source,[0 0]',x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),3e8);
x_isodop = fdoa.drawIsodop(x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),vdiff,1000,5);
plot(x_isodop(1,:),x_isodop(2,:),'-.','DisplayName','Line of Constant FDOA');
% text(.5,4.1,'FDOA Solution');

ylim([-5 5]);
xlim([-2 4]);
legend('Location','SouthWest');

utils.setPlotStyle(gca,{'clean','widescreen','tight'});
utils.exportPlot(fig4,[prefix '4']);


%% Figure 2.5 and 2.6
%  Example 2.1, generate data and 
addpath('examples');
figs = book2_ex2_1;

utils.exportPlot(figs(1),[prefix '5']);
utils.exportPlot(figs(2),[prefix '6a']);
utils.exportPlot(figs(3),[prefix '6b']);
utils.exportPlot(figs(4),[prefix '6c']);
utils.exportPlot(figs(5),[prefix '6d']);

%% Figure 2.7

fig = book2_ex2_2;

utils.exportPlot(fig, [prefix '7a']);

% Zoom in
xlim([2 3.5]*1e3);
ylim([2 3.5]*1e3);
utils.exportPlot(fig, [prefix '7b']);

%% Figure 2.8
%  Example 2.3 is lengthy to compute, skip this unless recalcluation of the
%  output is specifically requested.

if force_recalc
    fig = book2_ex2_3;

    utils.exportPlot(fig, [prefix '8a']);

    % Zoom in
    xlim([2 3.5]*1e3);
    ylim([2 3.5]*1e3);
    utils.exportPlot(fig, [prefix '8b']);
end

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;