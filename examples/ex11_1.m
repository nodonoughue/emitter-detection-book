function [fig_geo_a,fig_geo_b,fig_err] = ex11_1()
% [fig_geo_a,fig_geo_b,fig_err] = ex11_1()
%
% Executes Example 11.1 and generates two figures.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig_geo_a   figure handle for geographic layout
%   fig_geo_b   figure handle for geographic layout -- zoomed in on target
%   fig_err     figure handle for error as a function of iteration
%
% Nicholas O'Donoughue
% 1 July 2019

%  Set up TDOA Receiver system
%  Spacing of 1 km at 60 degree intervals around origin

% Define Sensor Positions
baseline = 10e3;
nSensors = 4;
thSensors = linspace(0,2*pi,nSensors) +pi/2; % add an extra sample, will be ignored
x_sensor = baseline* [ cos(thSensors(1:end-1));sin(thSensors(1:end-1))];

% Add one at the origin
x_sensor = cat(2,x_sensor,zeros(2,1));

% Define Sensor Performance
timingError = 1e-7;
c = 3e8;
rngStdDev = timingError*c;
Ctoa = timingError^2*eye(nSensors);
Croa = rngStdDev^2*eye(nSensors);
%Ctdoa = timingError^2*(1+eye(nSensors-1));
%Crdoa = rngStdDev^2*(1+eye(nSensors-1));

% Initialize Transmitter Position
th = rand(1)*2*pi;
x_source = 5*baseline*[cos(th);sin(th)];

% Compute Range
% R = sqrt(sum(abs(x_source-x_sensor).^2,1)); % Range to each sensor [m]
dR = tdoa.measurement(x_sensor,x_source);

% Generate Range Measurements
nMC = 1e3;
nr = rngStdDev*randn(nSensors,nMC); % [m]
ndr = nr(1:end-1,:) - nr(end,:); % Generate differential range msmnt noise
rho = dR + ndr; % Noisy sample vector

% Check covariance of rho
% mu_rho = mean(rho,2);
% cov_rho = cov(rho.');
% 
%     fprintf('Desired C:\n');
%     Crdoa
%     fprintf('Actual C:\n');
%     cov_rho
%     keyboard;
% 
x_init = 2*baseline*[cos(th+pi/6);sin(th+pi/6)]; % km

% Find Isochrones
xiso1 = tdoa.drawIsochrone(x_sensor(:,end),x_sensor(:,1),dR(1),1000,5*baseline);
xiso2 = tdoa.drawIsochrone(x_sensor(:,end),x_sensor(:,2),dR(2),1000,5*baseline);
xiso3 = tdoa.drawIsochrone(x_sensor(:,end),x_sensor(:,3),dR(3),1000,5*baseline);

% Call the TDOA solution algorithms
numIters = 1000;
alpha=.3;
beta=.8;
epsilon = 100; % desired search resolution size
x_ls = zeros(2,numIters);
thisX_ls = zeros(2,numIters);
x_grad = zeros(2,numIters);
x_chanHo = zeros(2,1);

cov_ls = zeros(2,2,numIters);
cov_grad = zeros(2,2,numIters);
cov_chanHo = zeros(2,2);
cov_bf = zeros(2,2);

fprintf('Performing Monte Carlo simulation...\n');
for idx = 1:nMC
    if mod(idx,floor(nMC/100))==0
        fprintf('.');
    end
    
    % Compute solutions
    [~,x_ls_local] = tdoa.lsSoln(x_sensor,rho(:,idx),Croa,x_init,[],numIters,true,[]);
    num_iter_ls = size(x_ls_local,2);
    if num_iter_ls > numIters
        x_ls_local = x_ls_local(:,1:numIters);
        num_iter_ls = numIters;
    end
    thisX_ls(:,1:num_iter_ls) = x_ls_local;
    if num_iter_ls < numIters
        thisX_ls(:,num_iter_ls+1:numIters) = x_ls_local(:,end);
    end
    
    [~,thisX_grad] = tdoa.gdSoln(x_sensor,rho(:,idx),Croa,x_init,alpha,beta,[],numIters,true,[]);
     
    thisX_chanHo = tdoa.chanHoSoln(x_sensor,rho(:,idx),Croa);
    
    thisX_bf = tdoa.bfSoln(x_sensor,rho(:,idx),Croa,x_init,5*baseline,epsilon);
    
    % Preserve first iteration for plotting
    if idx==1
        x_ls = thisX_ls;
        x_grad = thisX_grad;
        x_chanHo = thisX_chanHo;
        x_bf = thisX_bf;
    end
    
    % Compute Error
    thisErr_ls = x_source-thisX_ls; % [m]
    thisErr_grad = x_source-thisX_grad; % [m]
    thisErr_chanHo = x_source-thisX_chanHo;  % [m]
    thisErr_bf = x_source - thisX_bf; % [m]
    
    % Compute CEP50
    cov_ls = cov_ls + reshape(thisErr_ls,2,1,numIters).*reshape(thisErr_ls,1,2,numIters)/nMC;
    cov_grad = cov_grad + reshape(thisErr_grad,2,1,numIters).*reshape(thisErr_grad,1,2,numIters)/nMC;
    cov_chanHo = cov_chanHo + thisErr_chanHo*thisErr_chanHo'/nMC;
    cov_bf = cov_bf + thisErr_bf * thisErr_bf'/nMC;
end
fprintf('done.\n');

% Subsample iterative solutions for plotting
plotIndex = [1:10,20:20:100,200:200:numIters];
x_ls_plot = x_ls(:,plotIndex,1);
x_grad_plot = x_grad(:,plotIndex,1);



% Compute the covariance matrix for each iteration
cep50_ls = zeros(1,numIters);
cep50_grad = zeros(1,numIters);

for ii=1:numIters    
    cep50_ls(ii) = utils.computeCEP50(cov_ls(:,:,ii))/1e3; % [km]
    cep50_grad(ii) = utils.computeCEP50(cov_grad(:,:,ii))/1e3; % [km]
end
cep50_chanHo = utils.computeCEP50(cov_chanHo)/1e3; % [km]
cep50_bf = utils.computeCEP50(cov_bf)/1e3; % [km]

% Compute CRLB on RMSE
err_crlb = tdoa.computeCRLB(x_sensor,x_source,Ctoa);
crlb_cep50 = utils.computeCEP50(err_crlb)/1e3; % [km]
crlb_ellipse = utils.drawErrorEllipse(x_source,err_crlb,100,90);

% Plot Geometry
numIterToPlot=min(20,numel(plotIndex));
fig_geo_a=figure;
plot(x_sensor(1,:)/1e3,x_sensor(2,:)/1e3,'ko','DisplayName','Sensors','LineWidth',1);
hold on;
plot(x_source(1,:)/1e3,x_source(2,:)/1e3,'k^','DisplayName','Transmitter','LineWidth',1);
plot(x_ls_plot(1,1:floor(numIterToPlot/2),1)/1e3,x_ls_plot(2,floor(1:numIterToPlot/2),1)/1e3,':x','DisplayName',sprintf('Least Squares (%d iterations)',plotIndex(floor(numIterToPlot/2))),'LineWidth',1,'MarkerSize',4);
plot(x_grad_plot(1,1:numIterToPlot,1)/1e3,x_grad_plot(2,1:numIterToPlot,1)/1e3,':+','DisplayName',sprintf('Grad Descent (%d iterations)',plotIndex(numIterToPlot)),'LineWidth',1,'MarkerSize',4);
set(gca,'ColorOrderIndex',3);

% Label Solutions
text(x_init(1)/1e3+2,x_init(2)/1e3,'Initial Guess','FontSize',12);
text(mean(x_ls_plot(1,1:2,1))/1e3+5,mean(x_ls_plot(2,1:2,1))/1e3,'Least Squares','FontSize',12);
text(x_grad_plot(1,2,1)/1e3-40,x_grad_plot(2,2,1)/1e3,'Gradient Descent','FontSize',12);
annotation(fig_geo_a,'arrow',[0.612847222222222 0.559027777777778],...
    [0.643518518518518 0.5]);
annotation(fig_geo_a,'arrow',[0.665798611111111 0.733506944444444],...
    [0.600308641975309 0.387345679012346]);

xlabel('Cross-range [km]');ylabel('Down-range [km]');
legend('Location','NorthWest');
wd = 7*baseline/1e3;
ht = wd*9/16;
xlim([-1 1]*wd);
ylim([-1 1]*ht+x_source(2)/1e3/2);
grid off;
utils.setPlotStyle(gca,{'widescreen','equal','tight'});

% Plot zoomed geometry
fig_geo_b=figure;
plot(x_source(1,:)/1e3,x_source(2,:)/1e3,'k^','DisplayName','Transmitter');
hold on;
plot(x_ls_plot(1,1:floor(numIterToPlot/2),1)/1e3,x_ls_plot(2,1:floor(numIterToPlot/2),1)/1e3,'--x','DisplayName','LS');
plot(x_grad_plot(1,1:numIterToPlot,1)/1e3,x_grad_plot(2,1:numIterToPlot,1)/1e3,'-.+','DisplayName','Grad Descent');
plot(x_chanHo(1,1)/1e3,x_chanHo(2,1)/1e3,'*','DisplayName','Chan-Ho');
plot(x_bf(1,1)/1e3,x_bf(2,1)/1e3,'o','DisplayName','BestFix');

plot(crlb_ellipse(1,:)/1e3,crlb_ellipse(2,:)/1e3,'k','LineWidth',.5,'DisplayName','90% Error Ellipse');

wd = 1.4*max(abs(crlb_ellipse(1,:)-x_source(1)));
ht=wd*7/8;
ylim([-1 1]*ht/1e3+x_source(2)/1e3);
xlim([-1 1]*wd/1e3+x_source(1)/1e3);
xlabel('Cross-range [km]');ylabel('Down-range [km]');
legend('Location','best');
grid off;

utils.setPlotStyle(gca,{'equal','widescreen','tight'});

% Plot Error
fig_err=figure;
plot(1:numIters,cep50_ls,'DisplayName','LS');
hold on;
plot(1:numIters,cep50_grad,'--','DisplayName','Gradient Descent');
plot(1:numIters,cep50_chanHo*ones(1,numIters),'-.','DisplayName','Chan-Ho');
plot(1:numIters,cep50_bf*ones(1,numIters),'-+','DisplayName','BestFix');
plot(1:numIters,crlb_cep50*ones(1,numIters),':','DisplayName','CRLB CEP_{50}');
xlabel('Iteration Number');
ylabel('Position Error [km]');
set(gca,'yscale','log');
set(gca,'xscale','log');
text(22,15,'Gradient Descent','FontSize',10);
text(4,6,'Least Square','FontSize',10);
text(1.5,1.8,'Chan-Ho','FontSize',10);
text(1.5,1.2,'CRLB $CEP_{50}$','FontSize',10);
ylim([1 20]);
utils.setPlotStyle(gca,{'widescreen','tight'});
