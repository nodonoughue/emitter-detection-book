% Draw Figures - Chapter 1
%
% This script generates all of the figures that appear in
% Chapter 1 of the textbook.
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
prefix = fullfile(dirNm,'fig1_');

% Initialize Plot Preference
doFullColor=true;
utils.initPlotSettings;

% Reset the random number generator, to ensure reproducability
rng('default') ; 


%% Figure 1a, LOB Intersection
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
isoFill1 = cat(2,xaoap1,fliplr(xaoam1),xaoap1(:,1));

lob2 = x_source - x2;
aoa2 = atan2(lob2(2),lob2(1));
xaoa2 = x2 + [0 cos(aoa2);0 sin(aoa2)]*5*r2;
xaoap2 = x2 + [0 cos(aoa2+epsang);0 sin(aoa2+epsang)]*5*r2;
xaoam2 = x2 + [0 cos(aoa2-epsang);0 sin(aoa2-epsang)]*5*r2;
isoFill2 = cat(2,xaoap2,fliplr(xaoam2),xaoap2(:,1));

% Draw Figure
fig2a = figure();hold on;

% Uncertainty Intervals
fill(lobFill0(1,:),lobFill0(2,:),.2,'FaceAlpha',.1,'EdgeColor','k','LineStyle','--','DisplayName','Uncertainty Interval');
h = fill(isoFill1(1,:),isoFill1(2,:),.2,'FaceAlpha',.1,'EdgeColor','k','LineStyle','--');
utils.excludeFromLegend(h);
h = fill(isoFill2(1,:),isoFill2(2,:),.2,'FaceAlpha',.1,'EdgeColor','k','LineStyle','--');
utils.excludeFromLegend(h);


% LOBs
plot(xaoa0(1,:),xaoa0(2,:),'k-','DisplayName','AOA Solution');
h=plot(xaoa1(1,:),xaoa1(2,:),'k-');
h1=plot(xaoa2(1,:),xaoa2(2,:),'k-');
utils.excludeFromLegend([h,h1]);

% Position Markers
plot([x0(1),x1(1),x2(1)],[x0(2),x1(2),x2(2)],'ko','DisplayName','Sensors');
plot(x_source(1),x_source(2),'k^','MarkerSize',8,'DisplayName','Transmitter');

% Position Labels
text(x0(1)+.05,x0(2)-.1,'$S_0$');
text(x1(1)+.05,x1(2)-.1,'$S_1$');
text(x2(1)+.05,x2(2)-.1,'$S_2$');

% Adjust Axes
legend('Location','SouthEast');
utils.setPlotStyle(gca,{'clean','equal','widescreen','tight'});
ylim([-.5 1.5]);
xlim([-1 2]);

utils.exportPlot(fig2a,[prefix '2a']);

%% Figure 1b, TDOA Isochrone Intersection

% Initialize Detector/Source Locations
x_sensor1 = [0; 0];
x_sensor2 = [.8; .2];
x_sensor3 = [1;1];
x_source = [.1; .9];

% Compute Ranges
r1 = utils.rng(x_sensor1,x_source);
r2 = utils.rng(x_sensor2,x_source);
r3 = utils.rng(x_sensor3,x_source);

% Uncertainty Error
epsrdoa = .1;

% Find Isochrones
xiso1 = tdoa.drawIsochrone(x_sensor1,x_sensor2,r2-r1,1000,5);
xisop1 = tdoa.drawIsochrone(x_sensor1, x_sensor2, r2-r1+epsrdoa, 1000, 5);
xisom1 = tdoa.drawIsochrone(x_sensor1, x_sensor2, r2-r1-epsrdoa, 1000, 5);
isoFill1 = cat(2,xisop1,fliplr(xisom1),xisop1(:,1));

xiso2 = tdoa.drawIsochrone(x_sensor2,x_sensor3,r3-r2,1000,5);
xisop2 = tdoa.drawIsochrone(x_sensor2, x_sensor3, r3-r2+epsrdoa, 1000, 5);
xisom2 = tdoa.drawIsochrone(x_sensor2, x_sensor3, r3-r2-epsrdoa, 1000, 5);
isoFill2 = cat(2,xisop2,fliplr(xisom2),xisop2(:,1));

% Draw Figure
fig2b = figure();hold on;

% Uncertainty Intervals
fill(isoFill1(1,:),isoFill1(2,:),.2,'FaceAlpha',.1,'EdgeColor','k','LineStyle','--','DisplayName','Uncertainty Interval');
h = fill(isoFill2(1,:),isoFill2(2,:),.2,'FaceAlpha',.1,'EdgeColor','k','LineStyle','--');
utils.excludeFromLegend(h);

% Isochrones
plot(xiso1(1,:),xiso1(2,:),'k-','DisplayName','Isochrone');
hiso2=plot(xiso2(1,:),xiso2(2,:),'k-');
utils.excludeFromLegend(hiso2);

% Isochrone Labels
text(mean([x_sensor1(1),x_sensor2(1)])+.1,mean([x_sensor1(2),x_sensor2(2)])-.2,'$S_{1,2}$ Solution');
text(mean([x_sensor2(1),x_sensor3(1)])+.4,mean([x_sensor2(2),x_sensor3(2)])+.1,'$S_{2,3}$ Solution');

% Position Markers
plot([x_sensor1(1),x_sensor2(1),x_sensor3(1)],[x_sensor1(2),x_sensor2(2),x_sensor3(2)],'ko','DisplayName','Sensors');
plot(x_source(1),x_source(2),'k^','MarkerSize',8,'DisplayName','Transmitter');

% Position Labels
text(x_sensor1(1)+.05,x_sensor1(2)-.1,'$S_1$');
text(x_sensor2(1)+.05,x_sensor2(2)-.1,'$S_2$');
text(x_sensor3(1)+.05,x_sensor3(2)-.1,'$S_3$');

% Adjust Axes
xlim([-2 3]);
ylim([-1 2]);
legend('Location','SouthEast');
utils.setPlotStyle(gca,{'clean','widescreen','tight'});
utils.exportPlot(fig2b,[prefix '2b']);


%% Figure 1c, FDOA IsoDoppler Contour Intersection

x_source = [0;2]; % Transmitter/source
x_sensor = [-1,0;
       0,0;
       2,0]';
v_sensor = [1,1;
      1,1;
      1,1]';

% Draw Geometry
fig2c=figure;
hold on;

% Solve isodoppler line S12
vdiff12 = utils.dopDiff(x_source,[0 0]',x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),3e8);
xy_isodop12 = fdoa.drawIsodop(x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),vdiff12,1000,5);

xy_isodop12p = fdoa.drawIsodop(x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),vdiff12*1.1,1000,5);
xy_isodop12m = fdoa.drawIsodop(x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),vdiff12*0.9,1000,5);
fill12=cat(2,xy_isodop12p,fliplr(xy_isodop12m),xy_isodop12p(:,1));

% Solve isodoppler line S23
vdiff23 = utils.dopDiff(x_source,[0 0]',x_sensor(:,2),v_sensor(:,2),x_sensor(:,3),v_sensor(:,3),3e8);
xy_isodop23 = fdoa.drawIsodop(x_sensor(:,2),v_sensor(:,2),x_sensor(:,3),v_sensor(:,3),vdiff23,1000,5);

xy_isodop23p = fdoa.drawIsodop(x_sensor(:,2),v_sensor(:,2),x_sensor(:,3),v_sensor(:,3),vdiff23*1.1,1000,5);
xy_isodop23m = fdoa.drawIsodop(x_sensor(:,2),v_sensor(:,2),x_sensor(:,3),v_sensor(:,3),vdiff23*0.9,1000,5);
fill23=cat(2,xy_isodop23p,fliplr(xy_isodop23m),xy_isodop23p(:,1));

fill(fill12(1,:),fill12(2,:),.2,'FaceAlpha',.1,'EdgeColor','k','LineStyle','--','DisplayName','UncertaintyInterval');
h = fill(fill23(1,:),fill23(2,:),.2,'FaceAlpha',.1,'EdgeColor','k','LineStyle','--');
utils.excludeFromLegend(h);

% FDOA Solutions
plot(xy_isodop12(1,:),xy_isodop12(2,:),'k-','DisplayName','Line of Constant FDOA');
text(-1.1,2.5,'$S_{12}$ Solution','FontSize',10);

hh=plot(xy_isodop23(1,:),xy_isodop23(2,:),'k-');
utils.excludeFromLegend(hh);
text(1.6,.85,'$S_{23}$ Solution','FontSize',10);

% Position Markers
plot(x_source(1),x_source(2),'k^','DisplayName','Transmitter');
plot(x_sensor(1,:),x_sensor(2,:),'ko','DisplayName','Sensors');

% Draw velocity arrows
utils.drawArrow(x_sensor(1,1)+[0 v_sensor(1,1)]/4,x_sensor(2,1)+[0 v_sensor(2,1)]/4);
utils.drawArrow(x_sensor(1,2)+[0 v_sensor(1,2)]/4,x_sensor(2,2)+[0 v_sensor(2,2)]/4);
utils.drawArrow(x_sensor(1,3)+[0 v_sensor(1,3)]/4,x_sensor(2,3)+[0 v_sensor(2,3)]/4);

% Text labels
text(x_sensor(1,1)-.2,x_sensor(2,1)-.2,'$S_1$','FontSize',10);
text(x_sensor(1,2)-.2,x_sensor(2,2)-.2,'$S_2$','FontSize',10);
text(x_sensor(1,3)-.2,x_sensor(2,3)-.2,'$S_3$','FontSize',10);


xlim([-2 3]);
ylim([-1 4]);
legend;

utils.setPlotStyle(gca,{'clean','widescreen','tight'});
utils.exportPlot(fig2c,[prefix '2c']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;