% Draw Figures - Chapter 11
%
% This script generates all of the figures that appear in
% Chapter 11 of the textbook.
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
prefix = fullfile(dirNm,'fig11_');

% Initialize Plot Preference
utils.initPlotSettings;

% Initializes colorSet - Mx3 RGB vector for successive plot lines
colors = get(0,'DefaultAxesColorOrder');

% Reset the random number generator, to ensure reproducability
rng('default'); 

% Add path to folder for examples from the textbook
addpath('examples');

%% Figure 1a, TOA Circles

% Initialize Locations
x_sensor1 = [0; 0];
x_sensor2 = [1.5; .7];
x_sensor3 = [-.1; .8];
x_source = [.5; .5];

% Compute Range Vectors
r1 = utils.rng(x_sensor1,x_source);
r2 = utils.rng(x_sensor2,x_source);
r3 = utils.rng(x_sensor3,x_source);

% Generate a unit circle
th = linspace(0,2*pi,1001);
xcircle = cos(th);
ycircle = sin(th);
xx = [xcircle;ycircle]; % 2 x 1001

% Scale the unit circle and offset with position bias
toa1 = bsxfun(@plus,r1*xx,x_sensor1);
toa2 = bsxfun(@plus,r2*xx,x_sensor2);
toa3 = bsxfun(@plus,r3*xx,x_sensor3);

% Draw the circles and position hash marks
fig1a = figure();hold on;
plot(toa1(1,:),toa1(2,:),'k:','DisplayName','Constant TOA');
hiso2=plot(toa2(1,:),toa2(2,:),'k:',...
           toa3(1,:),toa3(2,:),'k:');
utils.excludeFromLegend(hiso2);

% Overlay Text
text(x_sensor1(1)+.05,x_sensor1(2)-.1,'$S_1$');
text(x_sensor2(1)+.05,x_sensor2(2)-.1,'$S_2$');
text(x_sensor3(1)+.05,x_sensor3(2)+.1,'$S_3$');

% text(x_source(1)+.1,x_source(2)+.1,'$T_1$');

% Plot radials
hiso2 = plot(x_sensor1(1)+[0,r1*cos(5*pi/4)],x_sensor1(2)+[0,r1*sin(5*pi/4)],'k-','LineWidth',.5);
utils.excludeFromLegend(hiso2);
plot(x_sensor2(1)+[0,r2*cos(5*pi/4)],x_sensor2(2)+[0,r2*sin(5*pi/4)],'k-','LineWidth',.5);
utils.excludeFromLegend(hiso2);
plot(x_sensor3(1)+[0,r3*cos(3*pi/4)],x_sensor3(2)+[0,r3*sin(3*pi/4)],'k-','LineWidth',.5);
utils.excludeFromLegend(hiso2);

text(x_sensor1(1)+r1/2*cos(5*pi/4)+.05,x_sensor1(2)+r1/2*sin(5*pi/4)-.1,'$R_1$');
text(x_sensor2(1)+r2/2*cos(5*pi/4)+.05,x_sensor2(2)+r2/2*sin(5*pi/4)-.1,'$R_2$');
text(x_sensor3(1)+r3/2*cos(3*pi/4),x_sensor3(2)+r3/2*sin(3*pi/4),'$R_3$');

% Position Markers (with legend)
plot([x_sensor1(1),x_sensor2(1),x_sensor3(1)],[x_sensor1(2),x_sensor2(2),x_sensor3(2)],'ko','DisplayName','Sensors');
plot(x_source(1),x_source(2),'k^','DisplayName','Transmitter');%,'MarkerSize',12);

% Adjust Axes
xlim([-2 3]);
ylim([-1 2]);
legend('Location','SouthEast');
utils.setPlotStyle(gca,{'clean','equal','widescreen','tight'});
utils.exportPlot(fig1a,[prefix '1a']);

%% Figure 1b TDOA Circles

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
fig1b = figure();hold on;

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
utils.setPlotStyle(gca,{'clean','widescreen','tight'});
utils.exportPlot(fig1b,[prefix '1b']);

%% Figure 2, Isochrones Plots

x_sensor1 = [0;0];
x_sensor2 = [1;0];

% Define set of RangeDiffs
Niso = 15;
rdiff = linspace(-.9,.9,Niso);

% Initialize Isochrone Matrices
Npts = 1000;
xiso = zeros(Niso,2*Npts-1);
yiso = zeros(Niso,2*Npts-1);

% Compute isochrones
for ind=1:Niso
    iso = tdoa.drawIsochrone(x_sensor1,x_sensor2,rdiff(ind),Npts,3);
    xiso(ind,:) = iso(1,:);
    yiso(ind,:) = iso(2,:);
end

% Initialize Figure
fig2 = figure();hold on;

% Plot isocrhones
hiso2 = plot(xiso',yiso','k-');
utils.excludeFromLegend(hiso2);

% Plot Markers
plot([x_sensor1(1),x_sensor2(1)],[x_sensor1(2),x_sensor2(2)],'ko','DisplayName','Sensors');

% Label the sensors
text(x_sensor1(1)-.1,x_sensor1(2),'$S_1$');
text(x_sensor2(1)+.05,x_sensor2(2),'$S_2$');

% Label the regimes
text(.24,2.2,'$R_1 < R_2$');
text(.55,2.2,'$R_1 > R_2$');
annotation(fig2,'arrow',[0.5-.0113 0.5-.0885],...
    [0.887345679012346 0.887345679012346]);
annotation(fig2,'arrow',[0.54+.0113 0.54+.0885],...
    [0.887345679012346 0.887345679012346]);

% Adjust Display
xlim([-0.5 1.5]);
ylim([-2.5 2.5]);
utils.setPlotStyle(gca,{'clean','widescreen','tight'});
utils.exportPlot(fig2,[prefix '2']);

%% Figure 3, Plot of leading edge detection

% Noise signal
s2 = .1;
N = 1024;
n = sqrt(s2)*randn(1,N);

% Signal
Es = 10;
t0 = floor(.3*N); % Start position of chirp
T = floor(N/2); % Length of chirp
p = sqrt(Es)*chirp(1:T,0,T,.1,'linear',-90);
y = zeros(1,N);
y(t0+(0:T-1)) = p;

fig3=figure;
plot(1:N,n,'DisplayName','Noise');
hold on;
set(gca,'ColorOrderIndex',3);
plot(1:N,y,'DisplayName','Signal');

% Lines
eta = 10*s2;
tt = find(abs(y)>eta,1,'first');
plot([1 N],eta*[1 1],'k--','DisplayName','Threshold');
plot(tt*[1 1],[-1 1]*sqrt(Es),'k:','DisplayName','\tau_i');
%legend('Location','NorthWest');

text(10,1.2,'Threshold');
text(10,-1.2,'Noise');
text(830,2.5,'Signal');
text(300,3,'$\tau_i$');
xlim([1 N]);
utils.setPlotStyle(gca,{'clean','widescreen','tight'});
utils.exportPlot(fig3,[prefix '3']);

%% Figure 4, Illustration of Cross Correlation TDOA Processing

% Noise signal
s2 = 5;
N = 1024;
n1 = sqrt(s2/2)*(randn(1,N)+1i*randn(1,N));
n2 = sqrt(s2/2)*(randn(1,N)+1i*randn(1,N));

% Signal
Es = 100;
t0 = floor(.3*N); % Start position of chirp
T = floor(N/2); % Length of chirp
p = sqrt(Es/2)*chirp(1:T,0,T,.1,'linear',-90);
y = zeros(1,N);
y(t0+(0:T-1)) = p+p*exp(1i*pi/2);

% Sensors Received Signals
s1 = n1+y;
s2 = n1+circshift(y,100);

% Time Index
fs = 1e4;
t = (0:N-1)/fs;
tau = t-mean(t);

sx = fftshift(ifft(fft(s1).*conj(fft(s2))))/N;

fig4=figure;
hold on;
subplot(311);hold on;
plot(t*1e3,real(s1));
title('Sensor 1');
plot(36*[1 1],[-10 10],'k:');
text(32,-5,'$\tau_1$');
ylim([-10 10]);
xlim([0 (N-1)*1e3/fs]);
utils.setPlotStyle(gca,'notick');
subplot(312);hold on;
plot(t*1e3,real(s2));
title('Sensor 2');
plot(46*[1 1],[-10 10],'k:');
text(42,-5,'$\tau_2$');
xlim([0 (N-1)*1e3/fs]);
ylim([-10 10]);
utils.setPlotStyle(gca,'notick');
subplot(313);hold on;
plot(tau*1e3,abs(sx));
plot(-10*[1 1],[0 40],'k:');
%plot([0 0],[0 40],'k:');
text(-22,30,'$\tau_{1,2}$');
%text(5,20,'0');
xlabel('Time Difference of Arrival [ms]');
title('Cross Correlation');
xlim([-1 1]*(N-1)*1e3/fs);
utils.setPlotStyle(gca,'notick');

utils.exportPlot(fig4,[prefix '4']);

%% Figure 5, Plot of the TDOA variance for peak detection and cross-correlation

% Link-16 Pulse Parameters
%   Chip Duration = 200 ns (implies 5 MHz bandwidth)
%     Assume flat pass-band shape
%   6.4 microseconds per pulse

% Define Pulse Parameters
bw = [5e6,1e9];
bw_rms = bw/sqrt(3);
pulseLen = [6.4e-6, 10e-6];
snrdb = linspace(0,60,100);
snrLin = 10.^(snrdb/10);

% Compute the variances
var_peak = tdoa.peakDetectionError(snrdb);
var_xcorr = tdoa.crossCorrError(snrdb,bw(:),pulseLen(:),bw_rms(:));

% Draw the figure
fig5 = figure();

semilogy(snrdb,sqrt(var_peak)*1e6,'DisplayName','Peak Detection');hold on;
txt = text(20,2e5,'Peak Detection');
set(txt,'Rotation',-8);
semilogy(snrdb,sqrt(var_xcorr(1,:))*1e6,'DisplayName',sprintf('Cross-Correlation, BT=%0.1f',pulseLen(1)*bw(1)));
txt = text(20,10,'Cross-Correlation, BT=32');
set(txt,'Rotation',-8);
semilogy(snrdb,sqrt(var_xcorr(2,:))*1e6,'DisplayName',sprintf('Cross-Correlalation, BT=%0.1f',pulseLen(2)*bw(2)));
txt = text(20,.02,'Cross-Correlation, BT=10,000');
set(txt,'Rotation',-8);

% Label Axes
ylabel('$\sigma$ [$\mu$s]');
xlabel('SNR [dB]');
%legend('Location','NorthEast');

utils.setPlotStyle(gca,{'widescreen','tight'});

% Save Output
utils.exportPlot(fig5,[prefix '5']);

%% Figure 6a
% Plots of the TDOA CRLB for a 3 sensor scenario and 1 microsecond timing 
% error

% Define Sensor Positions
baseline = 10e3;
nSensors = 3;
thSensors = linspace(0,2*pi,nSensors+1) +pi/2; % add an extra sample, will be ignored
x_sensor = baseline* [ cos(thSensors(1:end-1));sin(thSensors(1:end-1))];

% Define Sensor Performance
timingError = 1e-7;
Ctoa = timingError^2*eye(nSensors); % utilities now resample cov matrix with ref_idx
% Ctdoa = timingError^2 * (1 + eye(nSensors-1));

% Define source positions
M = 501;
xvec = linspace(-100,100,M)*1e3;
yvec = linspace(-100,100,M)*1e3;
[xx,yy] = ndgrid(xvec,yvec);
x_source = [xx(:) yy(:)]';

% Compute CRLB
warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
crlb = tdoa.computeCRLB(x_sensor,x_source,Ctoa); % Ndim x Ndim x M^2
cep50 = reshape(utils.computeCEP50(crlb),[M,M]);
warning('on','MATLAB:nearlySingularMatrix'); % Reactivate the singular matrix warning

% Set up contours
contourLevels=[.1,1,5,10,50,100,1000];
contourLevelsLabel = [.1,1,5,10,50,100];

% Draw Figure
fig6a = figure();hold on;

%ax=subplot(2,1,1)
plot(x_sensor(1,:)/1e3,x_sensor(2,:)/1e3,'ko','DisplayName','Sensors','LineWidth',1);
[cp,hiso2]=contour(xx/1e3,yy/1e3,cep50/1e3,contourLevels,'LineColor','k');
clabel(cp,hiso2);%contourLevelsLabel,fmt='%1.0f',fontsize=10,rightside_up=True)
utils.excludeFromLegend(hiso2);
legend('Location','NorthEast');
grid off;
% Adjust the Display
xlabel('Cross-range [km]');ylabel('Down-range [km]');
utils.setPlotStyle(gca,{'equal','tight'});

utils.exportPlot(fig6a,[prefix '6a']);

%% Figure 6b, Impact of fourth sensor on CRLB

% Add a sensor at the origin
x_sensor1 = [x_sensor zeros(2,1)];
nSensors = size(x_sensor1,2);

% Adjust Sensor Performance Vector
timingError = 1e-7;
Ctoa = timingError^2*eye(nSensors);
% Ctdoa = timingError^2 * (1 + eye(nSensors-1));

crlb2 = tdoa.computeCRLB(x_sensor1,x_source,Ctoa); % Ndim x Ndim x M**2
cep50 = reshape(utils.computeCEP50(crlb2),[M,M]);

% Draw the figure
fig6b = figure();hold on;

plot(x_sensor1(1,:)/1e3,x_sensor1(2,:)/1e3,'ko','DisplayName','Sensors','LineWidth',1);
[cp,hiso2]=contour(xx/1e3,yy/1e3,cep50/1e3,contourLevels,'LineColor','k');
clabel(cp,hiso2);%contourLevelsLabel,fmt='%1.0f',fontsize=10,rightside_up=True)
utils.excludeFromLegend(hiso2);
legend('Location','NorthEast')
grid off;

% Adjust the Display
xlabel('Cross-range [km]');ylabel('Down-range [km]');
utils.setPlotStyle(gca,{'equal','tight'});
utils.exportPlot(fig6b,[prefix '6b']);

%% Figures 7-8, Example TDOA Calculation
%  Figure 7 is geometry
%  Figure 8 is estimate error per iteration
if force_recalc
    
[fig7a,fig7b,fig8] = ex11_1;
utils.exportPlot(fig7a,[prefix '7a']);
utils.exportPlot(fig7b,[prefix '7b']);
utils.exportPlot(fig8,[prefix '8']);

end

%% Figure 9, Plot of false isochrones

% Define positions
x_sensor1 = [0;0];
x_sensor2 = [1;0];
x_source1 = [-.5;1];
x_source2 = [1;2];

% Ranges
r11 = utils.rng(x_source1,x_sensor1);
r12 = utils.rng(x_source1,x_sensor2);
r21 = utils.rng(x_source2,x_sensor1);
r22 = utils.rng(x_source2,x_sensor2);

% Compute Isochrones
Npts = 1000;
iso1 = tdoa.drawIsochrone(x_sensor1,x_sensor2,r12-r11,Npts,3);
iso2 = tdoa.drawIsochrone(x_sensor1,x_sensor2,r22-r21,Npts,3);
isoOff = tdoa.drawIsochrone(x_sensor1,x_sensor2,r22-r11,Npts,3);

% Draw Figure
fig9=figure();hold on;

% Transmitter/Sensor Locations
plot([x_sensor1(1),x_sensor2(1)],[x_sensor1(2),x_sensor2(2)],'ko','DisplayName','Sensors');
plot([x_source1(1),x_source2(1)],[x_source1(2),x_source2(2)],'k^','MarkerSize',8,'DisplayName','Transmitters');

% Transmitter/Sensor Labels
text(x_sensor1(1)-.4,x_sensor1(2)-.1,'$S_1$')
text(x_sensor2(1)+.1,x_sensor2(2)-.1,'$S_2$')

text(x_source1(1)+.1,x_source1(2)+.1,'$T_a$')
text(x_source2(1)+.1,x_source2(2)+.2,'$T_b$')


% Isochrones
hiso2 = plot(iso1(1,:),iso1(2,:),'k--');
utils.excludeFromLegend(hiso2);
plot(iso2(1,:),iso1(2,:),'k--','DisplayName','True Isochrones');
plot(isoOff(1,:),isoOff(2,:),'k:','DisplayName','False Isochrone');

% Isochrone Labels
labelYLoc = -1.5;
[~,ind1] = min(abs(iso1(2,:)-labelYLoc));
text(iso1(1,ind1)+.2,labelYLoc,'$TDOA_a$')
[~,ind2] = min(abs(iso2(2,:)-labelYLoc));
text(iso2(1,ind2)+.1,labelYLoc,'$TDOA_b$')
[~,indOff] = min(abs(isoOff(2,:)-labelYLoc));
text(isoOff(1,indOff)+.2,labelYLoc,'$TDOA_{err}$')

% Adjust Plot
legend('Location','West');
utils.setPlotStyle(gca,{'clean','widescreen','tight'});
utils.exportPlot(fig9,[prefix '9']);

%% Figure 10, Plot of false location

% Transmitter/Sensor Locations
x_sensor1 = [0;0];
x_sensor2 = [.8;.2];
x_sensor3 = [1;1];

x_source1 = [1;1.5];
x_source2 = [-.5;1.5];

% Ranges
r11 = utils.rng(x_source1,x_sensor1);
r12 = utils.rng(x_source1,x_sensor2);
r22 = utils.rng(x_source2,x_sensor2);
r23 = utils.rng(x_source2,x_sensor3);

% False Isochrones
Npts = 1000;
iso1 = tdoa.drawIsochrone(x_sensor1,x_sensor2,r12-r11,Npts,3);
iso2 = tdoa.drawIsochrone(x_sensor2,x_sensor3,r23-r22,Npts,3);

% Find False Solution
RR = utils.rng(iso1,iso2);
[~,idx] = min(RR(:));
[row,col] = ind2sub(size(RR),idx);
xFalse = mean([iso1(:,row),iso2(:,col)],2);

% Draw Figure
fig10=figure();hold on;

% Emiter/Sensor Locations
plot([x_sensor1(1),x_sensor2(1),x_sensor3(1)],[x_sensor1(2),x_sensor2(2),x_sensor3(2)],'ko','DisplayName','Sensors');
plot([x_source1(1),x_source2(1)],[x_source1(2),x_source2(2)],'k^','MarkerSize',8,'DisplayName','Transmitters');

% Transmitter/Sensor Labels
text(x_sensor1(1)+.05,x_sensor1(2)+.1,'$S_1$')
text(x_sensor2(1)+.05,x_sensor2(2)+.1,'$S_2$')
text(x_sensor3(1)+.05,x_sensor3(2)+.1,'$S_3$')

text(x_source1(1)+.1,x_source1(2),'$T_a$')
text(x_source2(1)+.1,x_source2(2),'$T_b$')

% False Isochrones
h1=plot(iso1(1,:),iso1(2,:),'k--');
h2=plot(iso2(1,:),iso2(2,:),'k--');
utils.excludeFromLegend([h1,h2]);
% Isochrone Labels
labelXLoc=2;
[~,ind1] = min(abs(iso1(1,:)-labelXLoc));
[~,ind2] = min(abs(iso2(1,:)-labelXLoc));
text(labelXLoc,iso1(2,ind1),'$TDOA_{12,a}$');
text(labelXLoc,iso2(2,ind2)+.2,'$TDOA_{23,b}$');

% False Solution
plot(xFalse(1),xFalse(2),'k^','MarkerFaceColor','k','MarkerSize',8,'DisplayName','False TDOA Solution')

% Adjust Display
legend('Location','SouthWest');
ylim([-1.5,2.5]);
xlim([-2.5,3.5]);
utils.setPlotStyle(gca,{'clean','widescreen','tight'});
utils.exportPlot(fig10,[prefix '10']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;