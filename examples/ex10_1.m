function [fig_geo,fig_err] = ex10_1()
% [fig_geo,fig_err] = ex10_1()
%
% Executes Example 10.1 and generates two figures
%
% INPUTS
%   none
%
% OUTPUTS
%   fig_geo     figure handle for geographic layout
%   fig_err     figure handle for error as a function of iteration
%
% Nicholas O'Donoughue
% 1 July 2019

% Define sensor positions
x_sensor = 30*[-1 0;0 0;1 0]';
N=size(x_sensor,2);

% Define source position
x_source = [15,45]';

% Grab a noisy measurement
psi_act = triang.measurement(x_sensor,x_source);
C_psi = (2*pi/180)^2*eye(N);
psi = psi_act + sqrt(C_psi)*randn(N,1);
r = 0:(1.2*max(utils.rng(x_source,x_sensor)));

% Compute Ranges
r0 = utils.rng(x_sensor(:,1),x_source);
r1 = utils.rng(x_sensor(:,2),x_source);
r2 = utils.rng(x_sensor(:,3),x_source);

% Error Values
epsang = 3*sqrt(diag(C_psi));

% Find AOA 
lob0 = x_source - x_sensor(:,1);
xaoa0 = x_sensor(:,1) + [0 cos(psi(1));0 sin(psi(1))]*5*r1;
xaoap0 = x_sensor(:,1) + [0 cos(psi(1)+epsang(1));0 sin(psi(1)+epsang(1))]*5*r1;
xaoam0 = x_sensor(:,1) + [0 cos(psi(1)-epsang(1));0 sin(psi(1)-epsang(1))]*5*r1;
lobFill0 = cat(2,xaoap0,fliplr(xaoam0),xaoap0(:,1));

lob1 = x_source - x_sensor(:,2);
xaoa1 = x_sensor(:,2) + [0 cos(psi(2));0 sin(psi(2))]*5*r1;
xaoap1 = x_sensor(:,2) + [0 cos(psi(2)+epsang(2));0 sin(psi(2)+epsang(2))]*5*r1;
xaoam1 = x_sensor(:,2) + [0 cos(psi(2)-epsang(2));0 sin(psi(2)-epsang(2))]*5*r1;
lobFill1 = cat(2,xaoap1,fliplr(xaoam1),xaoap1(:,1));

lob2 = x_source - x_sensor(:,3);
xaoa2 = x_sensor(:,3) + [0 cos(psi(3));0 sin(psi(3))]*5*r2;
xaoap2 = x_sensor(:,3) + [0 cos(psi(3)+epsang(3));0 sin(psi(3)+epsang(3))]*5*r2;
xaoam2 = x_sensor(:,3) + [0 cos(psi(3)-epsang(3));0 sin(psi(3)-epsang(3))]*5*r2;
lobFill2 = cat(2,xaoap2,fliplr(xaoam2),xaoap2(:,1));

% % Geometric Solutions
% x_centroid = triang.centroid(x_sensor,psi);
% x_incenter = triang.angle_bisector(x_sensor,psi);

% Iterative Methods
epsilon = .01;
nMC = 1e3;
numIters=100;
alpha=.3;
beta=.6;
% C_psi = eye(N);
psi = cat(2,psi,psi_act + sqrt(C_psi)*randn(N,nMC-1)); % preserve prior estimate
                                                       % so that geometric
                                                       % and iterative
                                                       % solutions use same
                                                       % inputs for plot

x_ml = zeros(2,nMC);          
x_bf = zeros(2,nMC);
x_centroid = zeros(2,nMC);
x_incenter = zeros(2,nMC);
x_ls_full = zeros(2,numIters,nMC);
x_grad_full = zeros(2,numIters,nMC);
fprintf('Conducting MC trial for triangulation error...\n');
for idx = 1:nMC
    if mod(idx,floor(nMC/80))==0
        fprintf('.');
    end
    x_ml(:,idx) = triang.mlSoln(x_sensor,psi(:,idx),C_psi,[5;5],[50;50],.1);
    x_bf(:,idx) = triang.bfSoln(x_sensor,psi(:,idx),C_psi,[5;5],[50;50],.1);
    x_centroid(:,idx) = triang.centroid(x_sensor,psi(:,idx));
    x_incenter(:,idx) = triang.angle_bisector(x_sensor,psi(:,idx));
    [~,x_ls_full(:,:,idx)] = triang.lsSoln(x_sensor,psi(:,idx),C_psi,[5;5],[],numIters,true,[]);
    [~,x_grad_full(:,:,idx)] = triang.gdSoln(x_sensor,psi(:,idx),C_psi,[5;5],[],[],[],numIters,true,[]);
end
err_ml = x_source - x_ml;
err_bf = x_source - x_bf;
err_cnt = x_source - x_centroid;
err_inc = x_source - x_incenter;
err_ls = x_source - x_ls_full;
err_grad = x_source - x_grad_full;

fprintf('done.\n');

bias_ml = mean(err_ml,2);
bias_bf = mean(err_bf,2);
bias_cnt = mean(err_cnt,2);
bias_inc = mean(err_inc,2);
cov_ml = cov(err_ml')+bias_ml*bias_ml';
cov_bf = cov(err_bf')+bias_bf*bias_bf';
cov_cnt = cov(err_cnt')+bias_cnt*bias_cnt';
cov_inc = cov(err_inc')+bias_inc*bias_inc';
cep50_ml = utils.computeCEP50(cov_ml);
cep50_bf = utils.computeCEP50(cov_bf);
cep50_cnt = utils.computeCEP50(cov_cnt);
cep50_inc = utils.computeCEP50(cov_inc);

bias_ls = zeros(2,numIters);
bias_grad = zeros(2,numIters);
cov_ls = zeros(2,2,numIters);
cov_grad = zeros(2,2,numIters);
cep50_ls = zeros(1,numIters);
cep50_grad = zeros(1,numIters);
for ii=1:numIters
    bias_ls(:,ii) = mean(squeeze(err_ls(:,ii,:)),2);
    bias_grad(:,ii) = mean(squeeze(err_grad(:,ii,:)),2);
    
    cov_ls(:,:,ii) = cov(squeeze(err_ls(:,ii,:))') + bias_ls(:,ii)*bias_ls(:,ii)';
    cov_grad(:,:,ii) = cov(squeeze(err_grad(:,ii,:))') + bias_grad(:,ii)*bias_grad(:,ii)';
    
    cep50_ls(ii) = utils.computeCEP50(cov_ls(:,:,ii)); % [km]
    cep50_grad(ii) = utils.computeCEP50(cov_grad(:,:,ii)); % [km]
end


%% First subfigure
% Draw Figure
fig_geo = figure();hold on;

% Uncertainty Intervals
h = fill(lobFill0(1,:),lobFill0(2,:),'k','FaceAlpha',.1,'EdgeColor','none');
utils.excludeFromLegend(h);

h = fill(lobFill1(1,:),lobFill1(2,:),'k','FaceAlpha',.1,'EdgeColor','none');
utils.excludeFromLegend(h);

h = fill(lobFill2(1,:),lobFill2(2,:),'k','FaceAlpha',.1,'EdgeColor','none');
utils.excludeFromLegend(h);

% Position Markers
hh=plot(x_sensor(1,:),x_sensor(2,:),'ko','DisplayName','Sensors');
utils.excludeFromLegend(hh);

% Position Labels
text(x_sensor(1,1)+2,x_sensor(2,1)-.1,'$S_0$');
text(x_sensor(1,2)+2,x_sensor(2,2)-.1,'$S_1$');
text(x_sensor(1,3)+2,x_sensor(2,3)-.1,'$S_2$');

% Plot the points and lobs
plot(x_source(1),x_source(2),'k^','DisplayName','Transmitter');

% Geometric Solutions
plot(x_centroid(1),x_centroid(2),'k*','DisplayName','Centroid');
plot(x_incenter(1),x_incenter(2),'k+','DisplayName','Incenter');
plot(x_ml(1,1),x_ml(2,1),'kv','DisplayName','Maximum Likelihood');
plot(x_bf(1,1),x_bf(2,1),'ko','DisplayName','BestFix');

% Iterative Solutions
plot(5,5,'kx','DisplayName','Initial Estimate');
hh=plot(x_ls_full(1,:,1),x_ls_full(2,:,10),'k:','DisplayName','Least Squares');
utils.excludeFromLegend(hh);
hh=plot(x_grad_full(1,:,1),x_grad_full(2,:,10),'k--','DisplayName','Grad Descent');
utils.excludeFromLegend(hh);
text(15,10,'Least Squares','FontSize',10);
text(-16,10,'Grad Descent','FontSize',10);
xlabel('[km]');
ylabel('[km]');

% Compute CRLB and Error Ellipse
err_crlb = triang.computeCRLB(x_sensor,x_source,C_psi);
crlb_cep50 = utils.computeCEP50(err_crlb); % [km]
crlb_ellipse = utils.drawErrorEllipse(x_source,err_crlb,100,90);
h=plot(crlb_ellipse(1,:),crlb_ellipse(2,:),'LineWidth',.5,'DisplayName','90% Error Ellipse');
utils.excludeFromLegend(h);
text(-20,45,'90\% Error Ellipse','FontSize',10);
h=plot([1 11],[45 45],'k-','LineWidth',.5);
utils.excludeFromLegend(h);
xlim([-50 50]);
ylim([-10 70]);

% Annotation Arrows
annotation(fig_geo,'arrow',[0.635451388888888 0.709236111111111],...
    [0.313924768518518 0.425601851851852]);
annotation(fig_geo,'arrow',[0.559027777777778 0.506944444444444],...
    [0.278645833333333 0.420572916666667]);

grid off;
legend('Location','NorthWest');
utils.setPlotStyle(gca,{'tight'});


%% Second subfigure
fig_err=figure;
loglog(1:numIters,cep50_ls,':','DisplayName','Least-Squares');hold on;
text(1.2,4,'Least-Squares','FontSize',10);
plot(1:numIters,cep50_grad,'--','DisplayName','Gradient Descent');
text(2.5,15,'Gradient Descent','FontSize',10);
plot(1:numIters,cep50_ml*ones(1,numIters),'DisplayName','Max Likelihood');
plot(1:numIters,cep50_bf*ones(1,numIters),'DisplayName','BestFix');
plot(1:numIters,cep50_cnt*ones(1,numIters),'DisplayName','Centroid');
text(15,2.6,'Centroid','FontSize',10);
plot(1:numIters,cep50_inc*ones(1,numIters),'DisplayName','Incenter');
text(3,35,'Incenter','FontSize',10);
plot(1:numIters,crlb_cep50*ones(1,numIters),'k-.','DisplayName','CRLB');
text(1.2,1.8,'CRLB','FontSize',10);

xlabel('Iteration Number');
ylabel('$CEP_{50}$ [km]');
% legend('location','NorthEast');
xlim([1 150]);
ylim([1 50]);

utils.setPlotStyle(gca,{'widescreen','tight'});
