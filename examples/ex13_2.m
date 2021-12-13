function [fig_geo,fig_err] = ex13_2()
% [fig_geo,fig_err] = ex13_2()
%
% Executes Example 13.2 and generates two figures.
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

x_source = [30;30]*1e3; % Transmitter/source
x_init = [0; 10e3];

x_aoa = 10e3 * [-1; 1];
n_aoa = 1;
x_tf = 10e3* [-1 1; 0 0];
v_tf = 100 * [0 0; 1 1];
n_tf = 2;

% Generate Measurements
c=3e8;
f_0 = 1e9;
z = hybrid.measurement(x_aoa,x_tf,x_tf,v_tf,x_source); % 3N x 1
ang_err = .06; % rad
time_err = 100e-9; % sec
rng_err = c*time_err; % m
freq_err = 3; % Hz
rng_rate_err = freq_err * c/f_0; % m/s
C_aoa = ang_err^2 * eye(n_aoa);
C_roa = rng_err^2 * eye(n_tf);
C_rroa = rng_rate_err^2 * eye(n_tf);
C_full = blkdiag(C_aoa, C_roa, C_rroa);

% Build reference vectors
aoa_test = 1:n_aoa;
aoa_ref = nan*ones(1,n_aoa);
tfdoa_ref = n_tf;
[tfdoa_test_vec, tfdoa_ref_vec] = utils.parseReferenceSensor(tfdoa_ref, n_tf);
test_vec = [aoa_test, n_aoa + tfdoa_test_vec, n_aoa + n_tf+ tfdoa_test_vec];
ref_vec = [aoa_ref, n_aoa + tfdoa_ref_vec, n_aoa + n_tf + tfdoa_ref_vec];

% Resample covariance matrix for tdoa/fdoa ref pairs
C_tilde = utils.resampleCovMtx(C_full, test_vec, ref_vec);

% Error Vectors
warning('Monte Carlo trials were set to 1,000 for the textbook, but 100 appears to be sufficient for stable results.  Lowered to 100 for faster execution.');
n_MC = 1e2; 
L = chol(C_tilde,'lower'); % Use the Cholesky decomposition to generate the square root
noise_full = L*randn(n_aoa+2*(n_tf-1),n_MC);

% Set up optimization parameters
num_iters = 1000;
plot_progress = false;
alpha = .3;
beta = .8;

% Perform Monte Carlo Loop
cov_ls = zeros(2,2,num_iters);
cov_grad = zeros(2,2,num_iters);
n_div_ls = 0;
n_div_grad = 0;

fprintf('Performing Monte Carlo simulation...\n');

warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
for idx_mc = 1:n_MC
    fprintf('.');
    if mod(idx_mc,50)==0
        fprintf(' (%d/%d)\n',idx_mc,n_MC);
    end

    % Add noise to measurement
    zeta = z+noise_full(:,idx_mc);

    % Least Squares
    [~,x_full] = hybrid.lsSoln(x_aoa,x_tf,x_tf,v_tf,zeta,C_full,x_init,[],num_iters,true,plot_progress,tfdoa_ref,tfdoa_ref);
    if idx_mc==1
        x_ls = x_full;
    end
    if any(isnan(x_full(:)))
        n_div_ls = n_div_ls+1;
    else
        err = x_full - x_source;
        cov_ls = cov_ls + reshape(err,2,1,num_iters).*conj(reshape(err,1,2,num_iters))/n_MC;
    end

    [~,x_full] = hybrid.gdSoln(x_aoa,x_tf,x_tf,v_tf,zeta,C_full,x_init,alpha,beta,[],num_iters,true,plot_progress,tfdoa_ref,tfdoa_ref);
    if idx_mc==1
        x_grad = x_full;
    end
    if any(isnan(x_full(:)))
        n_div_grad = n_div_grad+1;
    else
        err = x_full - x_source;
        cov_grad = cov_grad + reshape(err,2,1,num_iters).*conj(reshape(err,1,2,num_iters))/n_MC;
    end
end
warning('on','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
fprintf('done.\n');

% Adjust for divergent trials
cov_ls = cov_ls * n_MC/(n_MC-n_div_ls);
cov_grad = cov_grad * n_MC/(n_MC-n_div_grad);

% Subsample iterative solutions for plotting
plotIndex = [1:10,20:20:100,200:200:num_iters];
x_ls_plot = x_ls(:,plotIndex,1);
x_grad_plot = x_grad(:,plotIndex,1);

% Compute the CEP50 for each iteration
cep50_ls = zeros(1,num_iters);
cep50_grad = zeros(1,num_iters);
for ii=1:num_iters
    cep50_ls(ii) = utils.computeCEP50(cov_ls(:,:,ii))/1e3; % [km]
    cep50_grad(ii) = utils.computeCEP50(cov_grad(:,:,ii))/1e3; % [km]
end


% Compute CRLB on RMSE
err_crlb = hybrid.computeCRLB(x_aoa,x_tf,x_tf,v_tf,x_source,C_full,tfdoa_ref,tfdoa_ref);
crlb_cep50 = utils.computeCEP50(err_crlb)/1e3; % [km]
crlb_ellipse = utils.drawErrorEllipse(x_source,err_crlb,100,90);

% Plot Geometry
numIterToPlot=min(20,numel(plotIndex));
fig_geo=figure;
hold on;
set(gca,'ydir','normal');
plot(x_aoa(1,:)/1e3,x_aoa(2,:)/1e3,'ks','DisplayName','AOA Sensors','LineWidth',1);
hold on;
plot(x_tf(1,:)/1e3,x_tf(2,:)/1e3,'ko','DisplayName','TDOA/FDOA Sensors','LineWidth',1);
plot(x_source(1,:)/1e3,x_source(2,:)/1e3,'k^','DisplayName','Emitter','LineWidth',1);
plot(x_ls_plot(1,1:floor(numIterToPlot/2),1)/1e3,x_ls_plot(2,floor(1:numIterToPlot/2),1)/1e3,':x','DisplayName',sprintf('Least Squares (%d iterations)',plotIndex(floor(numIterToPlot/2))),'LineWidth',1,'MarkerSize',4);
plot(x_grad_plot(1,1:numIterToPlot,1)/1e3,x_grad_plot(2,1:numIterToPlot,1)/1e3,'-.+','DisplayName',sprintf('Grad Descent (%d iterations)',plotIndex(numIterToPlot)),'LineWidth',1,'MarkerSize',4);
set(gca,'ColorOrderIndex',3);
plot(crlb_ellipse(1,:)/1e3,crlb_ellipse(2,:)/1e3,'k','LineWidth',.5,'DisplayName','90% Error Ellipse');

% Draw Velocity Arrows
for idx_sensor = 1:n_tf
    thisX = x_tf(:,idx_sensor);
    thisV = v_tf(:,idx_sensor)/norm(v_tf(:,idx_sensor));
    utils.drawArrow(thisX(1)/1e3+[0 thisV(1)]*3,thisX(2)/1e3+[0 thisV(2)]*3);
end

% Annotation
annotation(fig_geo,'arrow',[0.339409722222222 0.439236111111111],...
    [0.429012345679012 0.510802469135802]);
annotation(fig_geo,'arrow',[0.359375 0.501736111111111],...
    [0.364197530864197 0.364197530864197]);

text(-7,13,'Least Square','FontSize',10)
text(3,6.5,'Gradient Descent','FontSize',10);
text(-8,9,'Initial Guess','FontSize',10);

legend('Location','NorthWest');
wd = 30;
ht = wd*9/16;
xlabel('Cross-range [km]');ylabel('Down-range [km]');
xlim([-1 1]*wd+x_source(1)/1e3/2);
ylim([-1 1]*ht+x_source(2)/1e3/2);
grid off;
utils.setPlotStyle(gca,{'widescreen','equal','tight'});

% Plot Error
fig_err=figure;
% plot(1:num_iters,cep50_ml*ones(1,num_iters),'DisplayName','Maximum Likelihood');hold on;
plot(1:num_iters,cep50_ls,'DisplayName','Least Squares');
hold on;
plot(1:num_iters,cep50_grad,'--','DisplayName','Gradient Descent');
plot(1:num_iters,crlb_cep50*ones(1,num_iters),':','DisplayName','CRLB CEP_{50}');
xlabel('Iteration Number');
ylabel('Position Error [km]');
set(gca,'yscale','log');
set(gca,'xscale','log');
% legend('Location','NorthEast');
ylim([.5 50])
text(22,20,'Gradient Descent','FontSize',10);
text(4.5,2.5,'Least Square','FontSize',10);
text(1.5,1.15,'CRLB $CEP_{50}$','FontSize',10);
utils.setPlotStyle(gca,{'widescreen','tight'});
