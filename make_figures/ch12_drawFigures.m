% Draw Figures - Chapter 12
%
% This script generates all of the figures that appear in
% Chapter 12 of the textbook.
%
% Nicholas O'Donoughue
% 1 July 2019


% Clear Figures
close all;

% Flag to force re-execution of long scripts
if ~exist('force_recalc','var')
    force_recalc = false;
end

% Set up directory and filenames for figures
dirNm = fullfile(pwd,'figures');
if ~exist(dirNm,'dir')
    mkdir(dirNm);
end
prefix = fullfile(dirNm,'fig12_');

% Initialize Plot Preference
utils.initPlotSettings;

% Initializes colorSet - Mx3 RGB vector for successive plot lines
colors = get(0,'DefaultAxesColorOrder');

% Reset the random number generator, to ensure reproducability
rng('default'); 

% Add path to the folder for examples from the textbook
addpath('examples');

%% Figure 1, System Drawing

x_source = [0;2]; % Transmitter/source
x_sensor = [-1,0;
       0,0;
       2,0]';
v_sensor = [1,1;
      1,1;
      1,1]';


% Draw Geometry
fig1=figure;
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
xy_isodop12 = fdoa.drawIsodop(x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),vdiff12,1000,5);
plot(xy_isodop12(1,:),xy_isodop12(2,:),'k-.','DisplayName','Line of Constant FDOA');
text(-1,2.7,'$S_{12}$ Solution','FontSize',10);

% Draw isodoppler line S23
vdiff23 = utils.dopDiff(x_source,[0 0]',x_sensor(:,2),v_sensor(:,2),x_sensor(:,3),v_sensor(:,3),3e8);
xy_isodop23 = fdoa.drawIsodop(x_sensor(:,2),v_sensor(:,2),x_sensor(:,3),v_sensor(:,3),vdiff23,1000,5);
hh=plot(xy_isodop23(1,:),xy_isodop23(2,:),'k-.');
utils.excludeFromLegend(hh);
text(1.5,.85,'$S_{23}$ Solution','FontSize',10);
xlim([-2 3]);
ylim([-1 4]);
legend;

utils.setPlotStyle(gca,{'widescreen','tight','box only'});
utils.exportPlot(fig1,[prefix '1']);

%% Figure 2, IsoDoppler Lines

x_sensor0 = [-2,0]';
x_sensor1 = [2,0]';
xx = -5:.01:5;
[X,Y] = meshgrid(xx);
xt = [X(:), Y(:)]';

v_sensor0 = [1,0]';
v_sensor1 = [1,0]';
vt = [0,0]';

f = 1e9;
ddop = utils.dopDiff(xt,vt,x_sensor0,v_sensor0,x_sensor1,v_sensor1,f);
fig2a=figure;contour(X,Y,reshape(ddop,size(X)),20,'k')
% colorbar;
%colormap(gray);
hold on;
plot(x_sensor0(1),x_sensor0(2),'ko','MarkerSize',10);
plot(x_sensor1(1),x_sensor1(2),'ko','MarkerSize',10);

% Draw Velocity Arrows
utils.drawArrow(x_sensor0(1)+[0 v_sensor0(1)],x_sensor0(2)+[0 v_sensor0(2)]);
utils.drawArrow(x_sensor1(1)+[0 v_sensor1(1)],x_sensor1(2)+[0 v_sensor1(2)]);

% Annotation Text
text(x_sensor0(1)-1,x_sensor0(2),'$S_0$','FontSize',10);
text(x_sensor1(1)-1,x_sensor1(2),'$S_1$','FontSize',10);
text(3.5,0,'$f_1 = f_0$','FontSize',10);
text(-5,0,'$f_1 = f_0$','FontSize',10);
text(-.75,0,'$f_1 > f_0$','FontSize',10);

utils.setPlotStyle(gca,{'clean','equal'});
utils.exportPlot(fig2a,[prefix '2a']);

v_sensor0 = [0,1]';
v_sensor1 = [0,1]';
vt = [0,0]';

f = 1e9;
ddop = utils.dopDiff(xt,vt,x_sensor0,v_sensor0,x_sensor1,v_sensor1,f);
fig2b=figure;contour(X,Y,reshape(ddop,size(X)),11,'k')
%colorbar;
%colormap(gray);
hold on;
plot(x_sensor0(1),x_sensor0(2),'ko','MarkerSize',10);
plot(x_sensor1(1),x_sensor1(2),'ko','MarkerSize',10);

% Draw Velocity Arrows
utils.drawArrow(x_sensor0(1)+[0 v_sensor0(1)],x_sensor0(2)+[0 v_sensor0(2)]);
utils.drawArrow(x_sensor1(1)+[0 v_sensor1(1)],x_sensor1(2)+[0 v_sensor1(2)]);

% Annotation Text
text(x_sensor0(1),x_sensor0(2)-1,'$S_0$','FontSize',10);
text(x_sensor1(1),x_sensor1(2)-1,'$S_1$','FontSize',10);
text(3,3,'$f_1 > f_0$','FontSize',10);
text(-3,3,'$f_1 < f_0$','FontSize',10);
text(-3,-3,'$f_1 > f_0$','FontSize',10);
text(3,-3,'$f_1 < f_0$','FontSize',10);


utils.setPlotStyle(gca,{'clean','equal'});
utils.exportPlot(fig2b,[prefix '2b']);

%% Figure 3 -- FDOA Error

x_source = [0 0;-1 1;1 0]';
v_source =[1 0;0 1;0 1]';
x_sensor = [1 3]';
% f0 = 1e9;

[eps1,x_vec,y_vec] = fdoa.fdoaErr(x_source(:,[1 2]),v_source(:,[1 2]),eye(2),x_sensor,4,1001);
[eps2,x_vec2,y_vec2] = fdoa.fdoaErr(x_source(:,[1 3]),v_source(:,[1 3]),eye(2),x_sensor,4,1001);
[eps3,x_vec3,y_vec3] = fdoa.fdoaErr(x_source(:,[2 3]),v_source(:,[2 3]),eye(2),x_sensor,4,1001);
[eps4,x_vec4,y_vec4] = fdoa.fdoaErr(x_source,v_source,eye(3),x_sensor,4,1001);

fig3a=figure;
imagesc(x_vec,y_vec,10*log10(eps1'));
colormap(gray);caxis([-80 0]+max(10*log10(eps1(:))));
hold on;
plot(x_source(1,1),x_source(2,1),'ko','MarkerSize',10);
plot(x_source(1,2),x_source(2,2),'ko','MarkerSize',10);
plot(x_sensor(1,1),x_sensor(2,1),'k^','MarkerSize',10);
set(gca,'ydir','normal');
grid off;

% Draw Velocity Arrows
utils.drawArrow(x_source(1,1)+[0 v_source(1,1)]/2,x_source(2,1)+[0 v_source(2,1)]/2);
utils.drawArrow(x_source(1,2)+[0 v_source(1,2)]/2,x_source(2,2)+[0 v_source(2,2)]/2);

% Annotation Text
text(x_source(1,1)-.6,x_source(2,1),'$S_1$','FontSize',10);
text(x_source(1,2)-.5,x_source(2,2)+.25,'$S_2$','FontSize',10);
text(x_sensor(1)+.25,x_sensor(2)-.25,'Source','FontSize',10);

utils.setPlotStyle(gca,{'equal','clean','tight'});
utils.exportPlot(fig3a,[prefix '3a']);

fig3b=figure;
imagesc(x_vec2,y_vec2,10*log10(eps2'));
colormap(gray);caxis([-80 0]+max(10*log10(eps2(:))));
hold on;
plot(x_source(1,1),x_source(2,1),'ko','MarkerSize',10);
plot(x_source(1,3),x_source(2,3),'ko','MarkerSize',10);
plot(x_sensor(1,1),x_sensor(2,1),'k^','MarkerSize',10);
set(gca,'ydir','normal');
grid off;

% Draw Velocity Arrows
utils.drawArrow(x_source(1,1)+[0 v_source(1,1)]/2,x_source(2,1)+[0 v_source(2,1)]/2);
utils.drawArrow(x_source(1,3)+[0 v_source(1,3)]/2,x_source(2,3)+[0 v_source(2,3)]/2);

% Annotation Text
text(x_source(1,1)-.6,x_source(2,1),'$S_1$','FontSize',10);
text(x_source(1,3)+.25,x_source(2,3)-.25,'$S_3$','FontSize',10);
text(x_sensor(1)+.25,x_sensor(2)-.25,'Source','FontSize',10);

utils.setPlotStyle(gca,{'equal','clean','tight'});
utils.exportPlot(fig3b,[prefix '3b']);

fig3c=figure;
imagesc(x_vec3,y_vec3,10*log10(eps3'));
colormap(gray);caxis([-80 0]+max(10*log10(eps3(:))));
hold on;
plot(x_source(1,2),x_source(2,2),'ko','MarkerSize',10);
plot(x_source(1,3),x_source(2,3),'ko','MarkerSize',10);
plot(x_sensor(1,1),x_sensor(2,1),'k^','MarkerSize',10);
set(gca,'ydir','normal');
grid off;

% Draw Velocity Arrows
utils.drawArrow(x_source(1,2)+[0 v_source(1,2)]/2,x_source(2,2)+[0 v_source(2,2)]/2);
utils.drawArrow(x_source(1,3)+[0 v_source(1,3)]/2,x_source(2,3)+[0 v_source(2,3)]/2);

% Annotation Text
text(x_source(1,2)-.5,x_source(2,2)+.25,'$S_2$','FontSize',10);
text(x_source(1,3)+.25,x_source(2,3)-.25,'$S_3$','FontSize',10);
text(x_sensor(1)+.25,x_sensor(2)-.25,'Source','FontSize',10);

utils.setPlotStyle(gca,{'equal','clean','tight'});
utils.exportPlot(fig3c,[prefix '3c']);




fig3d=figure;
imagesc(x_vec4,y_vec4,10*log10(eps4'));
colormap(gray);caxis([-80 0]+max(10*log10(eps4(:))));
hold on;
plot(x_source(1,1),x_source(2,1),'ko','MarkerSize',10);
plot(x_source(1,2),x_source(2,2),'ko','MarkerSize',10);
plot(x_source(1,3),x_source(2,3),'ko','MarkerSize',10);
plot(x_sensor(1,1),x_sensor(2,1),'k^','MarkerSize',10);
set(gca,'ydir','normal');
grid off;

% Draw Velocity Arrows
utils.drawArrow(x_source(1,1)+[0 v_source(1,1)]/2,x_source(2,1)+[0 v_source(2,1)]/2);
utils.drawArrow(x_source(1,2)+[0 v_source(1,2)]/2,x_source(2,2)+[0 v_source(2,2)]/2);
utils.drawArrow(x_source(1,3)+[0 v_source(1,3)]/2,x_source(2,3)+[0 v_source(2,3)]/2);

% Annotation Text
text(x_source(1,1)-.6,x_source(2,1),'$S_1$','FontSize',10);
text(x_source(1,2)-.6,x_source(2,2)+.25,'$S_2$','FontSize',10);
text(x_source(1,3)+.25,x_source(2,3)-.25,'$S_3$','FontSize',10);
text(x_sensor(1)+.25,x_sensor(2)-.25,'Source','FontSize',10);

utils.setPlotStyle(gca,{'equal','clean','tight'});
utils.exportPlot(fig3d,[prefix '3d']);

%% Figure 4, Two-Ship Configuration

% Plot a notional two-ship (over time)
x_sensor = [0 0;.5 -1]';
x1 = x_sensor + [0 3]';
x2 = x1 + [0 3]';
v_sensor = [0 1;0 1]';
v1 = v_sensor;
v2 = v1;

% Transmitter Position
x_source = [3 5]';

% Plot the Sensor Positions
fig4=figure;hold on;
colorSet = get(gca,'ColorOrder');
set(gca,'ColorOrderIndex',4);
plot(x_sensor(1,:),x_sensor(2,:),'o');
text(x_sensor(1,1)-.1,x_sensor(2,1)-.75,'$t_0$','FontSize',12,'Color',colorSet(4,:));
hh=utils.drawArrow(x_sensor(1,1)+[0 v_sensor(1,1)],x_sensor(2,1)+[0 v_sensor(2,1)]);
hh.Color = colorSet(4,:);
hh=utils.drawArrow(x_sensor(1,2)+[0 v_sensor(1,2)],x_sensor(2,2)+[0 v_sensor(2,2)]);
hh.Color = colorSet(4,:);

set(gca,'ColorOrderIndex',3);
plot(x1(1,:),x1(2,:),'o');
text(x1(1,1)-.1,x1(2,1)-.75,'$t_1$','FontSize',12,'Color',colorSet(3,:));
hh=utils.drawArrow(x1(1,1)+[0 v1(1,1)],x1(2,1)+[0 v1(2,1)]);
hh.Color = colorSet(3,:);
hh=utils.drawArrow(x1(1,2)+[0 v1(1,2)],x1(2,2)+[0 v1(2,2)]);
hh.Color = colorSet(3,:);

set(gca,'ColorOrderIndex',2);
plot(x2(1,:),x2(2,:),'o');
text(x2(1,1)-.1,x2(2,1)-.75,'$t_2$','FontSize',12,'Color',colorSet(2,:));
hh=utils.drawArrow(x2(1,1)+[0 v2(1,1)],x2(2,1)+[0 v2(2,1)]);
hh.Color = colorSet(2,:);
hh=utils.drawArrow(x2(1,2)+[0 v2(1,2)],x2(2,2)+[0 v2(2,2)]);
hh.Color = colorSet(2,:);

% set(gca,'ColorOrderIndex',1);
% plot(x3(1,:),x3(2,:),'o');
% hh=utils.drawArrow(x3(1,1)+[0 v3(1,1)],x3(2,1)+[0 v3(2,1)]);
% hh.Color = colorSet(1,:);
% hh=utils.drawArrow(x3(1,2)+[0 v3(1,2)],x3(2,2)+[0 v3(2,2)]);
% hh.Color = colorSet(1,:);

% Draw Isochrones

% Draw isodoppler line at time 0
vdiff = utils.dopDiff(x_source,[0 0]',x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),3e8);
xy_isodop = fdoa.drawIsodop(x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),vdiff,1000,20);
set(gca,'ColorOrderIndex',4);
hh=plot(xy_isodop(1,:),xy_isodop(2,:),'-.');
utils.excludeFromLegend(hh);

% Draw isodoppler line at time 1
vdiff = utils.dopDiff(x_source,[0 0]',x1(:,1),v1(:,1),x1(:,2),v1(:,2),3e8);
xy_isodop = fdoa.drawIsodop(x1(:,1),v1(:,1),x1(:,2),v1(:,2),vdiff,1000,20);
set(gca,'ColorOrderIndex',3);
hh=plot(xy_isodop(1,:),xy_isodop(2,:),'-.');
utils.excludeFromLegend(hh);

% Draw isodoppler line at time 2
vdiff = utils.dopDiff(x_source,[0 0]',x2(:,1),v2(:,1),x2(:,2),v2(:,2),3e8);
xy_isodop = fdoa.drawIsodop(x2(:,1),v2(:,1),x2(:,2),v2(:,2),vdiff,1000,20);
set(gca,'ColorOrderIndex',2);
hh=plot(xy_isodop(1,:),xy_isodop(2,:),'-.');
utils.excludeFromLegend(hh);

% % Draw isodoppler line at time 3
% vdiff = utils.dopDiff(xs,[0 0]',x3(:,1),v3(:,1),x3(:,2),v3(:,2),3e8);
% xy_isodop = fdoa.drawIsodop(x3(:,1),v3(:,1),x3(:,2),v3(:,2),vdiff,1000,20);
% set(gca,'ColorOrderIndex',1);
% hh=plot(xy_isodop(1,:),xy_isodop(2,:),'-.');
% utils.excludeFromLegend(hh);

plot(x_source(1),x_source(2),'k^');
text(x_source(1)+.1,x_source(2)-.75,'Transmitter','FontSize',12);

ylim([-5 15]);
xlim([-5 5]);
utils.setPlotStyle(gca,{'widescreen','tight','clean'});
utils.exportPlot(fig4,[prefix '4']);

%% Figure 5, Freq Estimation Error
snr_db = -10:.5:30;
snr_lin = 10.^(snr_db/10);

Ts = 1e-3;              % 1ms us pulse
Br = [1e3,1e6];         % 10 kHz receiver
T_samp = 1./(2*Br);     % Nyquist sampling rate
BT  = Ts.*Br;           % Time-Bandwidth product
N = Ts./T_samp;         % Nyquist sampling rate is twice the bandwidth

delta_f = @(Ts) 1./Ts;
sigma_f = fdoa.freqEstCRLB(T_samp(:),N(:),snr_db);
sigma_fd = fdoa.freqDiffEstCRLB(Ts,Br(:),snr_db);

fig5=figure;
hold on;
%plot(snr_db,delta_f(Ts)*ones(size(snr_db)),'k-.','DisplayName',sprintf('Nyquist, T_s=%d ms',Ts*1e3));
%plot(snr_db,sigma_f(N,snr_lin),'k--','DisplayName',sprintf('CRLB (freq), N=%d',N));
% plot(snr_db,sigma_f(10*N,snr_lin),'k--','DisplayName','CRLB (freq), 10 \mu s pulse');
for idx_bt = 1:numel(BT)
    set(gca,'ColorOrderIndex',2*(idx_bt-1)+1);
    plot(snr_db,2*sigma_f(idx_bt,:),'-.','LineWidth',2-idx_bt/2,'DisplayName',sprintf('\\sigma_{f}, BT=%d',BT(idx_bt)));
    set(gca,'ColorOrderIndex',2*(idx_bt-1)+1);
    plot(snr_db,sigma_fd(idx_bt,:),'--','LineWidth',2-idx_bt/2,'DisplayName',sprintf('\\sigma_{fd}, BT=%d',BT(idx_bt)));
end
set(gca,'yscale','log')
legend('Location','SouthWest')
xlabel('SNR [dB]');
ylabel('$\sigma_f$ [Hz]');

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig5,[prefix '5']);

%% Figure 6
% Plots of the FDOA CRLB for a 3 sensor scenario and 1 microsecond timing 
% error

% Define Sensor Positions
baseline = 10e3;
std_vel = 100;
nSensors = 3;
thSensors = linspace(0,2*pi,nSensors+1) +pi/2; % add an extra sample, will be ignored
x_sensor = baseline* [ cos(thSensors(1:end-1));sin(thSensors(1:end-1))];
v_sensor = 100 * [cos(thSensors(1:end-1));sin(thSensors(1:end-1))];

% Define Sensor Performance
freqError = 10; % 1 Hz resolution
c = 3e8;
f0 = 1e9;
rngRateStdDev = freqError*c/f0;
%C = eye(nSensors-1)+1; % covariance matrix structure
%Cfoa = freqError^2*ones(nSensors,1);
% Cfdoa = freqError^2*C;
Crroa = rngRateStdDev^2*eye(nSensors);
%Crrdoa = rngRateStdDev^2*C;

% Define source positions
M = 501;
xvec = linspace(-100,100,M)*1e3;
yvec = linspace(-100,100,M)*1e3;
[xx,yy] = ndgrid(xvec,yvec);
x_source = [xx(:) yy(:)]';

% Compute CRLB
warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
crlb = fdoa.computeCRLB(x_sensor,v_sensor,x_source,Crroa); % Ndim x Ndim x M^2
cep50 = reshape(utils.computeCEP50(crlb),[M,M]);
warning('on','MATLAB:nearlySingularMatrix'); % Reactivate the singular matrix warning

% Set up contours
contourLevels=[.1,1,5,10,50,100];
% contourLevelsLabel = [.1,1,5,10,50,100];

% Draw Figure
fig6a = figure();hold on;

%ax=subplot(2,1,1)
plot(x_sensor(1,:)/1e3,x_sensor(2,:)/1e3,'ko','DisplayName','Sensors','LineWidth',1);
[~,hiso2]=contour(xx/1e3,yy/1e3,cep50/1e3,contourLevels,'LineColor','k');
%clabel(cp,hiso2,contourLevelsLabel,fmt='%1.0f',fontsize=10,rightside_up=True)
utils.excludeFromLegend(hiso2);
legend('Location','NorthEast');
grid off;
% Adjust the Display
xlabel('Cross-range [km]');ylabel('Down-range [km]');
utils.setPlotStyle(gca,{'equal','tight'});

utils.exportPlot(fig6a,[prefix '6a']);

% Repeat with +x velocity
v_sensor = 100 * [ones(1,nSensors);zeros(1,nSensors)];

warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
crlb = fdoa.computeCRLB(x_sensor,v_sensor,x_source,Crroa); % Ndim x Ndim x M^2
cep50 = reshape(utils.computeCEP50(crlb),[M,M]);
warning('on','MATLAB:nearlySingularMatrix'); % Reactivate the singular matrix warning

% Draw Figure 6b
fig6b = figure();hold on;

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

utils.exportPlot(fig6b,[prefix '6b']);

%% Figure 6c,d, Impact of fourth sensor on CRLB

% Add a sensor at the origin
x1 = [x_sensor zeros(2,1)];
v_sensor = 100 * [cos(thSensors(1:end-1));sin(thSensors(1:end-1))];
v1 = cat(2,v_sensor,[1;0]);

nSensors = size(x1,2);

% Define Sensor Performance
freqError = 10; % 1 Hz resolution
c = 3e8;
f0 = 1e9;
rngRateStdDev = freqError*c/f0;
%C = eye(nSensors-1)+1; % covariance matrix structure
Cfoa = freqError^2*eye(nSensors,1);
%Cfdoa = freqError^2*C;
%Crrdoa = rngRateStdDev^2*C;
Crroa = rngRateStdDev^2*eye(nSensors);

% Define source positions
M = 501;
xvec = linspace(-100,100,M)*1e3;
yvec = linspace(-100,100,M)*1e3;
[xx,yy] = ndgrid(xvec,yvec);
x_source = [xx(:) yy(:)]';

% Compute CRLB
warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
crlb = fdoa.computeCRLB(x1,v1,x_source,Crroa); % Ndim x Ndim x M^2
cep50 = reshape(utils.computeCEP50(crlb),[M,M]);
warning('on','MATLAB:nearlySingularMatrix'); % Reactivate the singular matrix warning

% Set up contours
contourLevels=[.1,1,5,10,50,100];
contourLevelsLabel = [.1,1,5,10,50,100];

% Draw Figure
fig6c = figure();hold on;

%ax=subplot(2,1,1)
plot(x1(1,:)/1e3,x1(2,:)/1e3,'ko','DisplayName','Sensors','LineWidth',1);
[cp,hiso2]=contour(xx/1e3,yy/1e3,cep50/1e3,contourLevels,'LineColor','k');
clabel(cp,hiso2);%contourLevelsLabel,fmt='%1.0f',fontsize=10,rightside_up=True)
utils.excludeFromLegend(hiso2);
legend('Location','NorthEast');
grid off;
% Adjust the Display
xlabel('Cross-range [km]');ylabel('Down-range [km]');
utils.setPlotStyle(gca,{'equal','tight'});

utils.exportPlot(fig6c,[prefix '6c']);

% Repeat with +x velocity
v1 = 100 * [ones(1,nSensors);zeros(1,nSensors)];

warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
crlb = fdoa.computeCRLB(x1,v1,x_source,Crroa); % Ndim x Ndim x M^2
cep50 = reshape(utils.computeCEP50(crlb),[M,M]);
warning('on','MATLAB:nearlySingularMatrix'); % Reactivate the singular matrix warning

% Draw Figure
fig6d = figure();hold on;

%ax=subplot(2,1,1)
plot(x1(1,:)/1e3,x1(2,:)/1e3,'ko','DisplayName','Sensors','LineWidth',1);
[cp,hiso2]=contour(xx/1e3,yy/1e3,cep50/1e3,contourLevels,'LineColor','k');
clabel(cp,hiso2);%contourLevelsLabel,fmt='%1.0f',fontsize=10,rightside_up=True)
utils.excludeFromLegend(hiso2);
legend('Location','NorthEast');
grid off;
% Adjust the Display
xlabel('Cross-range [km]');ylabel('Down-range [km]');

utils.setPlotStyle(gca,{'equal','tight'});
utils.exportPlot(fig6d,[prefix '6d']);

%% Figures 7-8, Example FDOA Calculation
%  Figure 7 is geometry
%  Figure 8 is estimate error per iteration
if force_recalc
    
[fig7a,fig7b,fig8] = ex12_1;
utils.exportPlot(fig7a,[prefix '7a']);
utils.exportPlot(fig7b,[prefix '7b']);
utils.exportPlot(fig8,[prefix '8']);

end
%% Cleanup

% Restore plot settings
utils.resetPlotSettings;