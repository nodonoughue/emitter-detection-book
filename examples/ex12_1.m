function [fig_geo_a,fig_geo_b,fig_err] = ex12_1()
% [fig_geo_a, fig_geo_b, fig_err] = ex12_1()
%
% Executes Example 12.1 and generates three figures
%
% INPUTS
%   none
%
% OUTPUTS
%   fig_geo_a   figure handle for geographic layout
%   fig_geo_b   figure handle for geographic layout -- zomed in on target
%   fig_err     figure handle for error as a function of iteration
%
% Nicholas O'Donoughue
% 1 July 2019

% Set the random number generator to a predictable value, so that
% annotations are stable
rng(0824637);

%  Set up FDOA Receiver system
%  Spacing of 10 km at 60 degree intervals around origin

% Define Sensor Positions
baseline = 10e3;    % km
std_velocity = 100; % m/s
nSensors = 4;
thSensors = linspace(0,2*pi,nSensors) +pi/2; % add an extra sample, will be ignored
x_sensor = baseline * [ cos(thSensors(1:end-1));sin(thSensors(1:end-1))];
v_sensor = std_velocity * [ cos(thSensors(1:end-1));sin(thSensors(1:end-1))];

% Add one at the origin
x_sensor = cat(2,x_sensor,zeros(2,1));
v_sensor = cat(2,v_sensor,[std_velocity;0]); % Sensor at origin moving vertically

% Define Sensor Performance
freqError = 3; % 1 Hz resolution
c = 3e8;
f0 = 1e9;
rngRateStdDev = freqError*c/f0;
%C = eye(nSensors)+1; % covariance matrix structure
%Crrdoa = rngRateStdDev^2*C;
Crroa = rngRateStdDev^2*eye(nSensors);

% Initialize Transmitter Position
th = rand(1)*2*pi;
x_source = 5*baseline*[cos(th);sin(th)];

% Generate noise free measurement
r_dot_true = fdoa.measurement(x_sensor,v_sensor,x_source);

% Generate noisy measurements
nMC = 1e3;
nr = rngRateStdDev*randn(nSensors,nMC); % [m]
ndr = nr(1:end-1,:) - nr(end,:); % Generate differential range msmnt noise
rho_dot = r_dot_true + ndr; % Noisy sample vector

x_init = baseline * x_source./abs(x_source);

% Call the FDOA solution algorithms
numIters = 1000;
alpha=.3;
beta=.8;
epsilon=100;
x_ls = zeros(2,numIters);
thisX_ls = zeros(2,numIters);
x_grad = zeros(2,numIters);

cov_ls = zeros(2,2,numIters);
cov_grad = zeros(2,2,numIters);
cov_bf = zeros(2,2);

fprintf('Performing Monte Carlo simulation...\n');

warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
for idx = 1:nMC
    if mod(idx,floor(nMC/100))==0
        fprintf('.');
    end
    
    % Compute solutions
    [~,x_ls_local] = fdoa.lsSoln(x_sensor,v_sensor,rho_dot(:,idx),Crroa,x_init,[],numIters,true,[]);
    num_iter_ls = size(x_ls_local,2);
    if num_iter_ls > numIters
        x_ls_local = x_ls_local(:,1:numIters);
        num_iter_ls = numIters;
    end
    thisX_ls(:,1:num_iter_ls) = x_ls_local;
    if num_iter_ls < numIters
        thisX_ls(:,num_iter_ls+1:numIters) = x_ls_local(:,end);
    end
    
    [~,thisX_grad] = fdoa.gdSoln(x_sensor,v_sensor,rho_dot(:,idx),Crroa,x_init,alpha,beta,[],numIters,true,[]);
    
    thisX_bf = fdoa.bfSoln(x_sensor,v_sensor,rho_dot(:,idx),Crroa,x_init,5*baseline,epsilon);
    
    % Preserve first iteration for plotting
    if idx==1
        x_ls = thisX_ls;
        x_grad = thisX_grad;
        x_bf = thisX_bf;
    end
    
    % Compute Error
    thisErr_ls = x_source-thisX_ls; % [m]
    thisErr_grad = x_source-thisX_grad; % [m]
    thisErr_bf = x_source-thisX_bf; % [m]
    
    % Compute CEP50
    cov_ls = cov_ls + reshape(thisErr_ls,2,1,numIters).*reshape(thisErr_ls,1,2,numIters)/nMC;
    cov_grad = cov_grad + reshape(thisErr_grad,2,1,numIters).*reshape(thisErr_grad,1,2,numIters)/nMC;
    cov_bf = cov_bf + thisErr_bf * thisErr_bf'/nMC;
end
warning('on','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
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
cep50_bf = utils.computeCEP50(cov_bf)/1e3; % [km]

% Compute CRLB on RMSE
err_crlb = fdoa.computeCRLB(x_sensor,v_sensor,x_source,Crroa);
crlb_cep50 = utils.computeCEP50(err_crlb)/1e3; % [km]
crlb_ellipse = utils.drawErrorEllipse(x_source,err_crlb,100,90);

% Plot Geometry
numIterToPlot=min(20,numel(plotIndex));
fig_geo_a=figure;
plot(x_sensor(1,:)/1e3,x_sensor(2,:)/1e3,'ko','DisplayName','Sensors','LineWidth',1);
hold on;
plot(x_source(1,:)/1e3,x_source(2,:)/1e3,'k^','DisplayName','Transmitter','LineWidth',1);
plot(x_ls_plot(1,1:floor(numIterToPlot/2),1)/1e3,x_ls_plot(2,floor(1:numIterToPlot/2),1)/1e3,':x','DisplayName',sprintf('Least Squares (%d iterations)',plotIndex(floor(numIterToPlot/2))),'LineWidth',1,'MarkerSize',4);
plot(x_grad_plot(1,1:numIterToPlot,1)/1e3,x_grad_plot(2,1:numIterToPlot,1)/1e3,'-.+','DisplayName',sprintf('Grad Descent (%d iterations)',plotIndex(numIterToPlot)),'LineWidth',1,'MarkerSize',4);
set(gca,'ColorOrderIndex',3);
plot(crlb_ellipse(1,:)/1e3,crlb_ellipse(2,:)/1e3,'k','LineWidth',.5,'DisplayName','90% Error Ellipse');

% Draw Velocity Arrows
for idx_sensor = 1:size(x_sensor,2)
    thisX = x_sensor(:,idx_sensor);
    thisV = v_sensor(:,idx_sensor)/norm(v_sensor(:,idx_sensor));
    utils.drawArrow(thisX(1)/1e3+[0 thisV(1)]*5,thisX(2)/1e3+[0 thisV(2)]*5);
end

% Label Solutions
text(x_init(1)/1e3+2,x_init(2)/1e3,'Initial Estimate','FontSize',12);
text(x_ls_plot(1,2,1)/1e3+4,x_ls_plot(2,2,1)/1e3+5,'Least Squares','FontSize',12);
text(x_grad_plot(1,2,1)/1e3-40,x_grad_plot(2,2,1)/1e3+5,'Gradient Descent','FontSize',12);
annotation(fig_geo_a,'arrow',[0.584201388888889 0.607638888888889],...
    [0.643518518518518 0.479938271604938]);
annotation(fig_geo_a,'arrow',[0.542534722222222 0.488715277777778],...
    [0.683641975308642 0.578703703703704]);
xlabel('Cross-range [km]');ylabel('Down-range [km]');
legend('Location','SouthWest');
wd = 7*baseline/1e3;
ht = wd*9/16;
xlim([-1 1]*wd+x_source(1)/1e3/2);
ylim([-1 1]*ht+x_source(2)/1e3/2);
grid off;

utils.setPlotStyle(gca,{'widescreen','equal','tight'});

% Plot zoomed geometry
fig_geo_b=figure;
plot(x_source(1,:)/1e3,x_source(2,:)/1e3,'k^','DisplayName','Transmitter');
hold on;
plot(x_ls_plot(1,1:floor(numIterToPlot/2),1)/1e3,x_ls_plot(2,1:floor(numIterToPlot/2),1)/1e3,':x','DisplayName','LS');
plot(x_grad_plot(1,1:numIterToPlot,1)/1e3,x_grad_plot(2,1:numIterToPlot,1)/1e3,'-.+','DisplayName','Grad Descent');
plot(x_bf(1,1)/1e3,x_bf(2,1)/1e3,'o','DisplayName','BestFix');
plot(crlb_ellipse(1,:)/1e3,crlb_ellipse(2,:)/1e3,'k','LineWidth',.5,'DisplayName','90% Error Ellipse');

wd = 8*max(abs(crlb_ellipse(1,:)-x_source(1)));
ht=wd*9/16;
ylim([-1 1]*ht/1e3+x_source(2)/1e3);
xlabel('Cross-range [km]');ylabel('Down-range [km]');
legend('Location','best');
grid off;

utils.setPlotStyle(gca,{'widescreen','equal','tight'});

% Plot Error
fig_err=figure;
plot(1:numIters,cep50_ls,'DisplayName','LS');
hold on;
plot(1:numIters,cep50_grad,'--','DisplayName','Gradient Descent');
plot(1:numIters,cep50_bf*ones(1,numIters),'-.','DisplayName','BestFix');
plot(1:numIters,crlb_cep50*ones(1,numIters),':','DisplayName','CRLB CEP_{50}');
xlabel('Iteration Number');
ylabel('Position Error [km]');
set(gca,'yscale','log');
set(gca,'xscale','log');
text(13,12,'Gradient Descent','FontSize',10);
text(3,4,'Least Square','FontSize',10);
text(1.05,1.75,'CRLB $CEP_{50}$','FontSize',10);
xlim([1,numIters]);
utils.setPlotStyle(gca,{'widescreen','tight'});
