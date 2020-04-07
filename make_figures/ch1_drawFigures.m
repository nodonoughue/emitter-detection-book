% Draw Figures - Chapter 1
%
% This script generates all of the figures that appear in
% Chapter 1 of the textbook.
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
prefix = fullfile(dirNm,'fig1_');

% Initialize Plot Preference
utils.initPlotSettings;

% Reset the random number generator, to ensure reproducability
rng('default') ; 


%% Figure 1, Detection Threshold
N0 = 0; % dB
M = 512; % number of points
n = sqrt(10^(N0/10)/2)* (randn(1,M)+1i*randn(1,M));

% Threshold
PFA = 1e-12;
T = sqrt(-log(PFA));

% Manually spike one noise sample
mspike = randi(M);
n(mspike) = T+1;

% Target Positions
m1 = fix(M/3);
a1 = T + 3;
m2 = fix(5*M/8);
a2 = T - 1;

% Target signal - use the autocorrelation of a window of length N to
% generate the lobing structure
t = 10*pi*linspace(-1,1,M);
p = sinc(t)/sqrt(sum(sinc(t).^2));
s = zeros(1,M);
s = s + circshift(p,m1)*10^(a1/10);
s = s + circshift(p,m2)*10^(a2/10);
s = ifft(fft(s,M).*conj(fft(p,M)),M);

% Plot Noise and Threshold
fig1=figure;hold on;
plot(1:M,10*log10(abs(s)),'DisplayName','Signals','linewidth',2);
plot(1:M,10*log10(abs(n)),'DisplayName','Noise','linewidth',.5);
plot([1 M],T*[1 1],'--','DisplayName','Threshold');
ylim([-5,T+5]);

text(m1+(M/50),a1,'Detection','FontSize',12);
text(m2,a2+.5,'Missed Detection','FontSize',12);
text(mspike-10,T+4,'False Alarm','FontSize',12);

xlim([1 M]);
legend('Location','NorthEast');

utils.setPlotStyle(gca,{'clean','widescreen','tight'});
utils.exportPlot(fig1,[prefix '1']);

%% Figure 2, AOA Geometry
% Compute an AOA slice from sensor 1

% Initialize Detector/Source Locations
x1 = [0; 0];
xs = [.1; .9];

% Compute Ranges
r1 = utils.rng(x1,xs);

% Error Values
epsang = 5*pi/180;

% Find AOA 
lob = xs - x1;
aoa1 = atan2(lob(2),lob(1));
xaoa1 = x1 + [0 cos(aoa1);0 sin(aoa1)]*5*r1;
xaoap1 = x1 + [0 cos(aoa1+epsang);0 sin(aoa1+epsang)]*5*r1;
xaoam1 = x1 + [0 cos(aoa1-epsang);0 sin(aoa1-epsang)]*5*r1;
lobFill1 = cat(2,xaoap1,fliplr(xaoam1),xaoap1(:,1));

% Draw Figure
fig2 = figure();hold on;

% LOBs
plot(xaoa1(1,:),xaoa1(2,:),'k-','DisplayName','AOA Solution');

% Uncertainty Intervals
plot(xaoap1(1,:),xaoap1(2,:)','k--','DisplayName','Uncertainty Interval');
h=plot(xaoam1(1,:),xaoam1(2,:),'k--');
utils.excludeFromLegend(h);
h = fill(lobFill1(1,:),lobFill1(2,:),'k','FaceAlpha',.1,'EdgeColor','none');
utils.excludeFromLegend(h);

% Position Markers
plot(x1(1),x1(2),'ko','DisplayName','Sensor');
plot(xs(1),xs(2),'k^','MarkerSize',8,'DisplayName','Transmitter');

% Adjust Axes
legend('Location','SouthEast');
utils.setPlotStyle(gca,{'clean','equal','widescreen','tight'});
ylim([-.5 1.5]);
utils.exportPlot(fig2,[prefix '2']);

%% Figure 3, Geolocation Geometry
% Compute an isochrone between sensors 1 and 2, then draw an AOA slice from
% sensor 3

% Initialize Detector/Source Locations
x1 = [0; 0];
x2 = [1;1];
xs = [.1; .9];

% Compute Ranges
r1 = utils.rng(x1,xs);
r2 = utils.rng(x2,xs);

% Error Values
epsang = 5*pi/180;

% Find AOA 
lob1 = xs - x1;
aoa1 = atan2(lob1(2),lob1(1));
xaoa1 = x1 + [0 cos(aoa1);0 sin(aoa1)]*5*r2;
xaoap1 = x1 + [0 cos(aoa1+epsang);0 sin(aoa1+epsang)]*5*r2;
xaoam1 = x1 + [0 cos(aoa1-epsang);0 sin(aoa1-epsang)]*5*r2;
lobFill1 = cat(2,xaoap1,fliplr(xaoam1),xaoap1(:,1));

lob2 = xs - x2;
aoa2 = atan2(lob2(2),lob2(1));
xaoa2 = x2 + [0 cos(aoa2);0 sin(aoa2)]*5*r2;
xaoap2 = x2 + [0 cos(aoa2+epsang);0 sin(aoa2+epsang)]*5*r2;
xaoam2 = x2 + [0 cos(aoa2-epsang);0 sin(aoa2-epsang)]*5*r2;
lobFill2 = cat(2,xaoap2,fliplr(xaoam2),xaoap2(:,1));

% Draw Figure
fig3 = figure();hold on;

% LOBs
plot(xaoa1(1,:),xaoa1(2,:),'k-','DisplayName','AOA Solution');
h=plot(xaoa2(1,:),xaoa2(2,:),'k-');
utils.excludeFromLegend(h);

% Uncertainty Intervals
plot(xaoap1(1,:),xaoap1(2,:)','k--','DisplayName','Uncertainty Interval');
h=plot(xaoam1(1,:),xaoam1(2,:),'k--');
utils.excludeFromLegend(h);
h = fill(lobFill1(1,:),lobFill1(2,:),'k','FaceAlpha',.1,'EdgeColor','none');
utils.excludeFromLegend(h);

h=plot(xaoap2(1,:),xaoap2(2,:)','k--',xaoam2(1,:),xaoam2(2,:),'k--');
utils.excludeFromLegend(h);
h = fill(lobFill2(1,:),lobFill2(2,:),'k','FaceAlpha',.1,'EdgeColor','none');
utils.excludeFromLegend(h);

% Position Markers
plot([x1(1),x2(1)],[x1(2),x2(2)],'ko','DisplayName','Sensors');
plot(xs(1),xs(2),'k^','MarkerSize',8,'DisplayName','Transmitter');

% Position Labels
text(x1(1)+.05,x1(2)-.1,'$S_1$');
text(x2(1)+.05,x2(2)-.1,'$S_2$');

% Adjust Axes
legend('Location','SouthEast');
utils.setPlotStyle(gca,{'clean','equal','widescreen','tight'});
ylim([-.5 1.5]);
xlim([-1 2]);

utils.exportPlot(fig3,[prefix '3']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;