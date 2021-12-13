% Draw Figures - Chapter 13
%
% This script generates all of the figures that appear in
% Chapter 13 of the textbook.
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
prefix = fullfile(dirNm,'fig13_');

% Initialize Plot Preference
utils.initPlotSettings;

% Initializes colorSet - Mx3 RGB vector for successive plot lines
colors = get(0,'DefaultAxesColorOrder');

% Reset the random number generator, to ensure reproducability
rng('default'); 

% Add folder for textbook examples
addpath('examples');

c=utils.constants.c;

%% Figure 1, System Drawing

x_source = [2;4]; % Transmitter/source

x_sensor = [0,0; 3,0]'; % Receivers
v_sensor = [1,0;1,0]';

fig1=figure;
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
psi = atan2(lob(2,:),lob(1,:));

xaoa1 = x_sensor(:,1) + [0 cos(psi(1));0 sin(psi(1))]*5*r(1);
xaoap1 = x_sensor(:,1) + [0 cos(psi(1)+epsang);0 sin(psi(1)+epsang)]*5*r(1);
xaoam1 = x_sensor(:,1) + [0 cos(psi(1)-epsang);0 sin(psi(1)-epsang)]*5*r(1);
lobFill1 = cat(2,xaoap1,fliplr(xaoam1),xaoap1(:,1));

xaoa2 = x_sensor(:,2) + [0 cos(psi(2));0 sin(psi(2))]*5*r(2);
xaoap2 = x_sensor(:,2) + [0 cos(psi(2)+epsang);0 sin(psi(2)+epsang)]*5*r(2);
xaoam2 = x_sensor(:,2) + [0 cos(psi(2)-epsang);0 sin(psi(2)-epsang)]*5*r(2);
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
xy_isodop = fdoa.drawIsodop(x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),vdiff,1000,5);
plot(xy_isodop(1,:),xy_isodop(2,:),'-.','DisplayName','Line of Constant FDOA');
% text(.5,4.1,'FDOA Solution');

ylim([-5 5]);
xlim([-2 4]);
legend('Location','SouthWest');

% Export Plot
utils.setPlotStyle(gca,{'widescreen','tight','clean'});
utils.exportPlot(fig1,[prefix '1']);

%% Figure 2 - ML Plot

x_sensor = [0 0;3 0]'*10e3;
v_sensor =[1 0;1 0]'*100;

x_source = [2 4]'*10e3;
f_0 = 10e9;

plot_dynamic_range = [-100 0];
contour_levels = -100:20:0;
% Generate the noise-free measurement
zeta = hybrid.measurement(x_sensor,x_sensor,x_sensor,v_sensor,x_source);

% Generate grid
grid_vec = -50e3:500:50e3;
[XX,YY] = meshgrid(grid_vec);
x_grid = [XX(:)+x_source(1) YY(:)+x_source(2)]';

% --- Figure 2a --- Baseline
% Generate noise covariance
n_sensor = size(x_sensor,2);
ang_err = .2; % rad
time_err = 1000e-9; % sec
rng_err = c*time_err; % m
freq_err = 100; % Hz
rng_rate_err = freq_err * c/f_0; % m/s

% Error Covariance Matrices
C_psi = (ang_err)^2 * eye(n_sensor);
C_roa = rng_err^2 * eye(n_sensor);
C_rroa = rng_rate_err^2 * eye(n_sensor);
C_full = blkdiag(C_psi,C_roa,C_rroa);

% Compute log likelihood
ella = hybrid.loglikelihood(x_sensor,x_sensor,x_sensor,v_sensor,zeta,C_full,x_grid);
ella = reshape(ella,size(XX));

% --- Figure 2b --- Better DF
% Generate noise covariance
ang_err = .02; % deg
time_err = 1000e-9; % sec
rng_err = c*time_err; % m
freq_err = 100; % Hz
rng_rate_err = freq_err * c/f_0; % m/s

% Error Covariance Matrices
C_psi = (ang_err)^2 * eye(n_sensor);
C_roa = rng_err^2 * eye(n_sensor);
C_rroa = rng_rate_err^2 * eye(n_sensor);
C_full = blkdiag(C_psi,C_roa,C_rroa);

% Compute log likelihood
ellb = hybrid.loglikelihood(x_sensor,x_sensor,x_sensor,v_sensor,zeta,C_full,x_grid);
ellb = reshape(ellb,size(XX));

% --- Figure 2c --- Better TDOA
% Generate noise covariance
ang_err = .2; % rad
time_err = 100e-9; % sec
rng_err = c*time_err; % m
freq_err = 100; % Hz
rng_rate_err = freq_err * c/f_0; % m/s

% Error Covariance Matrices
C_psi = (ang_err)^2 * eye(n_sensor);
C_roa = rng_err^2 * eye(n_sensor);
C_rroa = rng_rate_err^2 * eye(n_sensor);
C_full = blkdiag(C_psi,C_roa,C_rroa);

% Compute log likelihood
ellc = hybrid.loglikelihood(x_sensor,x_sensor,x_sensor,v_sensor,zeta,C_full,x_grid);
ellc = reshape(ellc,size(XX));

% --- Figure 2d --- Better FDOA
% Generate noise covariance
n_sensors = size(x_sensor,2);
ang_err = .2; % rad
time_err = 1000e-9; % sec
rng_err = c*time_err; % m
freq_err = 10; % Hz
rng_rate_err = freq_err * c/f_0; % m/s

% Error Covariance Matrices
C_psi = (ang_err)^2 * eye(n_sensor);
C_roa = rng_err^2 * eye(n_sensor);
C_rroa = rng_rate_err^2 * eye(n_sensor);
C_full = blkdiag(C_psi,C_roa,C_rroa);

% Compute log likelihood
elld = hybrid.loglikelihood(x_sensor,x_sensor,x_sensor,v_sensor,zeta,C_full,x_grid);
elld = reshape(elld,size(XX));

% --- Plot it all
fig2a=figure;
% imagesc(x_grid(1,:),x_grid(2,:),ella);
% colormap(flipud(gray));
% caxis(max(ella(:))+plot_dynamic_range);
[cc,hh] = contour((XX+x_source(1))/1e3,...
                  (YY+x_source(2))/1e3,ella,contour_levels,'k');
clabel(cc,hh,'Color','k');
utils.excludeFromLegend(hh);
hold on;
plot(x_sensor(1,:)/1e3,x_sensor(2,:)/1e3,'ko','MarkerSize',8,'DisplayName','Sensors');
plot(x_source(1,1)/1e3,x_source(2,1)/1e3,'k^','MarkerSize',8,'DisplayName','Source');
set(gca,'ydir','normal');
grid off;

% Draw Velocity Arrows
utils.drawArrow(x_sensor(1,1)/1e3+[0 v_sensor(1,1)]/50,x_sensor(2,1)/1e3+[0 v_sensor(2,1)]/50);
utils.drawArrow(x_sensor(1,2)/1e3+[0 v_sensor(1,2)]/50,x_sensor(2,2)/1e3+[0 v_sensor(2,2)]/50);
xlabel('Cross-Range [km]');
ylabel('Range [km]');
legend('Location','NorthWest');

utils.setPlotStyle(gca,{'widescreen','equal','tight'});
utils.exportPlot(fig2a,[prefix '2a']);

fig2b=figure;
% imagesc(x_grid(1,:),x_grid(2,:),ellb);
% colormap(flipud(gray));
% caxis(max(ellb(:))+plot_dynamic_range);
[cc,hh] = contour((XX+x_source(1))/1e3,...
                  (YY+x_source(2))/1e3,ellb,contour_levels,'k');
% clabel(cc,hh,'Color','k');
utils.excludeFromLegend(hh);
hold on;
plot(x_sensor(1,:)/1e3,x_sensor(2,:)/1e3,'ko','MarkerSize',8,'DisplayName','Sensors');
plot(x_source(1,1)/1e3,x_source(2,1)/1e3,'k^','MarkerSize',8,'DisplayName','Source');
set(gca,'ydir','normal');
grid off;

% Draw Velocity Arrows
utils.drawArrow(x_sensor(1,1)/1e3+[0 v_sensor(1,1)]/50,x_sensor(2,1)/1e3+[0 v_sensor(2,1)]/50);
utils.drawArrow(x_sensor(1,2)/1e3+[0 v_sensor(1,2)]/50,x_sensor(2,2)/1e3+[0 v_sensor(2,2)]/50);
xlabel('Cross-Range [km]');
ylabel('Range [km]');
legend('Location','NorthWest');

utils.setPlotStyle(gca,{'widescreen','equal','tight'});
utils.exportPlot(fig2b,[prefix '2b']);

fig2c=figure;
% imagesc(xs(1)+grid_vec,xs(2)+grid_vec,ellc);
% colormap(flipud(gray));
% caxis(max(ellc(:))+plot_dynamic_range);
[cc,hh] = contour((XX+x_source(1))/1e3,...
                  (YY+x_source(2))/1e3,ellc,contour_levels,'k');
% clabel(cc,hh,'Color','k');
utils.excludeFromLegend(hh);
hold on;
plot(x_sensor(1,:)/1e3,x_sensor(2,:)/1e3,'ko','MarkerSize',8,'DisplayName','Sensors');
plot(x_source(1,1)/1e3,x_source(2,1)/1e3,'k^','MarkerSize',8,'DisplayName','Source');
set(gca,'ydir','normal');
grid off;

% Draw Velocity Arrows
utils.drawArrow(x_sensor(1,1)/1e3+[0 v_sensor(1,1)]/50,x_sensor(2,1)/1e3+[0 v_sensor(2,1)]/50);
utils.drawArrow(x_sensor(1,2)/1e3+[0 v_sensor(1,2)]/50,x_sensor(2,2)/1e3+[0 v_sensor(2,2)]/50);
xlabel('Cross-Range [km]');
ylabel('Range [km]');
legend('Location','NorthWest');

utils.setPlotStyle(gca,{'widescreen','equal','tight'});
utils.exportPlot(fig2c,[prefix '2c']);


fig2d=figure;
% imagesc(x_grid(1,:),x_grid(2,:),elld);
% colormap(flipud(gray));
% caxis(max(elld(:))+plot_dynamic_range);
[cc,hh] = contour((XX+x_source(1))/1e3,...
                  (YY+x_source(2))/1e3,elld,contour_levels,'k');
% clabel(cc,hh,'Color','k');
utils.excludeFromLegend(hh);
hold on;
plot(x_sensor(1,:)/1e3,x_sensor(2,:)/1e3,'ko','MarkerSize',8,'DisplayName','Sensors');
plot(x_source(1,1)/1e3,x_source(2,1)/1e3,'k^','MarkerSize',8,'DisplayName','Source');
set(gca,'ydir','normal');
grid off;

% Draw Velocity Arrows
utils.drawArrow(x_sensor(1,1)/1e3+[0 v_sensor(1,1)]/50,x_sensor(2,1)/1e3+[0 v_sensor(2,1)]/50);
utils.drawArrow(x_sensor(1,2)/1e3+[0 v_sensor(1,2)]/50,x_sensor(2,2)/1e3+[0 v_sensor(2,2)]/50);
xlabel('Cross-Range [km]');
ylabel('Range [km]');
legend('Location','NorthWest');

utils.setPlotStyle(gca,{'widescreen','equal','tight'});
utils.exportPlot(fig2d,[prefix '2d']);


%% Figures 3,4 - Homogeneous (3-Mode) Sensors
if force_recalc
    
[fig3,fig4] = ex13_1;

utils.exportPlot(fig3,[prefix '3']);
utils.exportPlot(fig4,[prefix '4']);
end

%% Figures 5,6 - Heterogeneous Sensors
if force_recalc
    
[fig5,fig6] = ex13_2;
utils.exportPlot(fig5,[prefix '5']);
utils.exportPlot(fig6,[prefix '6']);

end

%% Figure 7 -- CRLB
% Define Sensor Positions
baseline = 10e3;
std_vel = 100;
n_sensor = 2;
x_sensor = baseline*[-.5 0;.5 0]';
v_sensor = std_vel*[1 0;1 0]';

% Define Sensor Performance
c=3e8;
f_0=1e9;
ang_err = .06;
time_err = 1e-7; % 100 ns resolution
rng_err = c*time_err; % m
freq_err = 10; % Hz
rng_rate_err = freq_err * c/f_0; % m/s

% Error Covariance Matrices
C_psi = (ang_err)^2 * eye(n_sensor);
C_roa = rng_err^2 * eye(n_sensor);
C_rroa = rng_rate_err^2 * eye(n_sensor);
C_full = blkdiag(C_psi,C_roa,C_rroa);

% Define source positions
M = 501;
xvec = linspace(-100,100,M)*1e3;
yvec = linspace(-100,100,M)*1e3;
[xx,yy] = ndgrid(xvec,yvec);
x_source = [xx(:) yy(:)]';

% Compute CRLB
warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
crlb = hybrid.computeCRLB(x_sensor,x_sensor,x_sensor,v_sensor,x_source,C_full); % Ndim x Ndim x M^2
cep50 = reshape(utils.computeCEP50(crlb),[M,M]);
warning('on','MATLAB:nearlySingularMatrix'); % Reactivate the singular matrix warning

% Set up contours
contourLevels=[.1,1,5,10];
contourLevelsLabel = [.1,1,5,10];

% Draw Figure
fig7a = figure();hold on;

%ax=subplot(2,1,1)
plot(x_sensor(1,:)/1e3,x_sensor(2,:)/1e3,'ko','DisplayName','Sensors','LineWidth',1);
[cp,hiso2]=contour(xx/1e3,yy/1e3,cep50/1e3,contourLevels,'LineColor','k');
clabel(cp,hiso2);%contourLevelsLabel,fmt='%1.0f',fontsize=10,rightside_up=True)
utils.excludeFromLegend(hiso2);
% Draw Velocity Arrows
for idx_sensor = 1:size(x_sensor,2)
    thisX = x_sensor(:,idx_sensor);
    thisV = v_sensor(:,idx_sensor)/norm(v_sensor(:,idx_sensor));
    utils.drawArrow(thisX(1)/1e3+[0 thisV(1)]*5,thisX(2)/1e3+[0 thisV(2)]*5);
end

legend('Location','NorthEast');
grid off;
% Adjust the Display
xlabel('Cross-range [km]');ylabel('Down-range [km]');

utils.setPlotStyle(gca,{'equal','tight'});

utils.exportPlot(fig7a,[prefix '7a']);

% Repeat with +x velocity
v_sensor = std_vel*[0 1;0 1]';

warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
crlb = hybrid.computeCRLB(x_sensor,x_sensor,x_sensor,v_sensor,x_source,C_full); % Ndim x Ndim x M^2
cep50 = reshape(utils.computeCEP50(crlb),[M,M]);
warning('on','MATLAB:nearlySingularMatrix'); % Reactivate the singular matrix warning

% Draw Figure
fig7b = figure();hold on;

%ax=subplot(2,1,1)
plot(x_sensor(1,:)/1e3,x_sensor(2,:)/1e3,'ko','DisplayName','Sensors','LineWidth',1);
[cp,hiso2]=contour(xx/1e3,yy/1e3,cep50/1e3,contourLevels,'LineColor','k');
clabel(cp,hiso2);%contourLevelsLabel,fmt='%1.0f',fontsize=10,rightside_up=True)
utils.excludeFromLegend(hiso2);
% Draw Velocity Arrows
for idx_sensor = 1:size(x_sensor,2)
    thisX = x_sensor(:,idx_sensor);
    thisV = v_sensor(:,idx_sensor)/norm(v_sensor(:,idx_sensor));
    utils.drawArrow(thisX(1)/1e3+[0 thisV(1)]*5,thisX(2)/1e3+[0 thisV(2)]*5);
end
legend('Location','NorthEast');
grid off;
% Adjust the Display
xlabel('Cross-range [km]');ylabel('Down-range [km]');

utils.setPlotStyle(gca,{'equal','tight'});

utils.exportPlot(fig7b,[prefix '7b']);

%% Figure 8 -- TDOA/FDOA only
% Define Sensor Positions
baseline = 10e3;
std_vel = 100;
n_sensor = 2;
x_sensor = baseline*[-.5 0;.5 0]';
v_sensor = std_vel*[1 0;1 0]';

% Define Sensor Performance
c=3e8;
f_0=1e9;
time_err = 1e-7; % 100 ns resolution
freq_err = 10; % 1 Hz resolution
rng_err = c*time_err; % m
rng_rate_err = freq_err * c/f_0; % m/s

% Error Covariance Matrices
C_roa = rng_err^2 * eye(n_sensor);
C_rroa = rng_rate_err^2 * eye(n_sensor);
C_full = blkdiag(C_roa,C_rroa);

% Define source positions
M = 501;
xvec = linspace(-100,100,M)*1e3;
yvec = linspace(-100,100,M)*1e3;
[xx,yy] = ndgrid(xvec,yvec);
x_source = [xx(:) yy(:)]';

% Compute CRLB
warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
crlb = hybrid.computeCRLB([],x_sensor,x_sensor,v_sensor,x_source,C_full); % Ndim x Ndim x M^2
cep50 = reshape(utils.computeCEP50(crlb),[M,M]);
warning('on','MATLAB:nearlySingularMatrix'); % Reactivate the singular matrix warning

% Set up contours
contourLevels=[.1,1,5,10,50];
contourLevelsLabel = [.1,1,5,10,50];

% Draw Figure
fig8a = figure();hold on;

%ax=subplot(2,1,1)
plot(x_sensor(1,:)/1e3,x_sensor(2,:)/1e3,'ko','DisplayName','Sensors','LineWidth',1);
[cp,hiso2]=contour(xx/1e3,yy/1e3,cep50/1e3,contourLevels,'LineColor','k');
clabel(cp,hiso2);%contourLevelsLabel,fmt='%1.0f',fontsize=10,rightside_up=True)
utils.excludeFromLegend(hiso2);
% Draw Velocity Arrows
for idx_sensor = 1:size(x_sensor,2)
    thisX = x_sensor(:,idx_sensor);
    thisV = v_sensor(:,idx_sensor)/norm(v_sensor(:,idx_sensor));
    utils.drawArrow(thisX(1)/1e3+[0 thisV(1)]*5,thisX(2)/1e3+[0 thisV(2)]*5);
end

legend('Location','NorthEast');
grid off;
% Adjust the Display
xlabel('Cross-range [km]');ylabel('Down-range [km]');

utils.setPlotStyle(gca,{'equal','tight'});

utils.exportPlot(fig8a,[prefix '8a']);

% Repeat with +x velocity
v_sensor = std_vel*[0 1;0 1]';

warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
crlb = hybrid.computeCRLB([],x_sensor,x_sensor,v_sensor,x_source,C_full); % Ndim x Ndim x M^2
cep50 = reshape(utils.computeCEP50(crlb),[M,M]);
warning('on','MATLAB:nearlySingularMatrix'); % Reactivate the singular matrix warning

% Draw Figure
fig8b = figure();hold on;

%ax=subplot(2,1,1)
plot(x_sensor(1,:)/1e3,x_sensor(2,:)/1e3,'ko','DisplayName','Sensors','LineWidth',1);
[cp,hiso2]=contour(xx/1e3,yy/1e3,cep50/1e3,contourLevels,'LineColor','k');
clabel(cp,hiso2);%contourLevelsLabel,fmt='%1.0f',fontsize=10,rightside_up=True)
utils.excludeFromLegend(hiso2);
% Draw Velocity Arrows
for idx_sensor = 1:size(x_sensor,2)
    thisX = x_sensor(:,idx_sensor);
    thisV = v_sensor(:,idx_sensor)/norm(v_sensor(:,idx_sensor));
    utils.drawArrow(thisX(1)/1e3+[0 thisV(1)]*5,thisX(2)/1e3+[0 thisV(2)]*5);
end
legend('Location','NorthEast');
grid off;
% Adjust the Display
xlabel('Cross-range [km]');ylabel('Down-range [km]');

utils.setPlotStyle(gca,{'equal','tight'});

utils.exportPlot(fig8b,[prefix '8b']);

%% Figure 9 -- TDOA/AOA only
% Define Sensor Positions
baseline = 10e3;
std_vel = 100;
n_sensor = 2;
x_sensor = baseline*[-.5 0;.5 0]';

% Define Sensor Performance
c=3e8;
ang_err = .06;
time_err = 1e-7; % 100 ns resolution
rng_err = c*time_err; % m

% Error Covariance Matrices
C_psi = (ang_err)^2 * eye(n_sensor);
C_roa = rng_err^2 * eye(n_sensor);
C_full = blkdiag(C_psi,C_roa);

% Define source positions
M = 501;
xvec = linspace(-100,100,M)*1e3;
yvec = linspace(-100,100,M)*1e3;
[xx,yy] = ndgrid(xvec,yvec);
x_source = [xx(:) yy(:)]';

% Compute CRLB
warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
crlb = hybrid.computeCRLB(x_sensor,x_sensor,[],[],x_source,C_full); % Ndim x Ndim x M^2
cep50 = reshape(utils.computeCEP50(crlb),[M,M]);
warning('on','MATLAB:nearlySingularMatrix'); % Reactivate the singular matrix warning

% Set up contours
contourLevels=[.1,1,5,10,50];
contourLevelsLabel = [.1,1,5,10,50];

% Draw Figure
fig9 = figure();hold on;

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

utils.exportPlot(fig9,[prefix '9']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;