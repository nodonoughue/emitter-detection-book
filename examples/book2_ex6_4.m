function figs = book2_ex6_4()
% figs = book2_ex6_4()
%
% Executes Example 6.4 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   figs         array of figure handle
%
% Nicholas O'Donoughue
% 12 February 2022

%% Set up sensors
x_tdoa = [2, 0,  4, 0;
          2, 2, 0, 0]*1e3; % avg position (reported)
n_tdoa = size(x_tdoa,2);

%% Generate Measurements
x_tgt = [6; 3]*1e3;

alpha_tdoa = [10, 30, -20, 60]'; % TOA bias
z = tdoa.measurement(x_tdoa,x_tgt, n_tdoa, alpha_tdoa); % free of pos unc and bias
z_true = tdoa.measurement(x_tdoa,x_tgt, n_tdoa); % free of pos unc and bias

err_toa = 100e-9;
C_toa = err_toa^2*eye(n_tdoa);
C_roa = utils.constants.c^2*C_toa;
C_rdoa = utils.resampleCovMtx(C_roa,n_tdoa);
U= chol(C_rdoa,'lower');

noise = U*randn(size(U,2),1);
zeta = z + noise;
zeta_true = z_true + noise;

%% Compute Log Likelihood
x_vec_km = 0:.05:10;
[xx,yy] = meshgrid(x_vec_km);
x_grid = [xx(:),yy(:)]'*1e3; % meters

ell = tdoa.loglikelihood(x_tdoa,zeta,C_roa,x_grid,n_tdoa);
ell_true = tdoa.loglikelihood(x_tdoa,zeta_true,C_roa,x_grid,n_tdoa);

%% Plot Scenario
fig1 = figure;
ell_true_plot = reshape(ell_true,size(xx))-max(ell_true);
levels=[-1000,-100,-50,-20,-10,-5,0];
[C,h]=contourf(x_vec_km,x_vec_km,ell_true_plot,levels);
clabel(C,h);
utils.excludeFromLegend(h);
% imagesc(x_vec_km,x_vec_km,ell_true_plot); % alternative -- use imagesc to plot
hold on;set(gca,'ydir','normal');
scatter(x_tdoa(1,:)/1e3,x_tdoa(2,:)/1e3,'^','filled','DisplayName','Sensors');
scatter(x_tgt(1)/1e3,x_tgt(2)/1e3,'o','filled','DisplayName','Target');
colorbar;
colormap(utils.viridis);
caxis([-100,0]);
grid on
legend('Location','NorthWest');

utils.setPlotStyle(gca,{'tight'});

fig2 = figure;
ell_plot = reshape(ell,size(xx))-max(ell);
[C,h]=contourf(x_vec_km,x_vec_km,ell_plot,levels);
clabel(C,h);
utils.excludeFromLegend(h);
hold on;set(gca,'ydir','normal');
scatter(x_tdoa(1,:)/1e3,x_tdoa(2,:)/1e3,'^','filled','DisplayName','Sensors');
scatter(x_tgt(1)/1e3,x_tgt(2)/1e3,'o','filled','DisplayName','Target');
colorbar;
colormap(utils.viridis);
caxis([-100,0]);
grid on;
legend('Location','Northwest');

utils.setPlotStyle(gca,{'tight'});

%% ML Solver (baseline)
x_ctr = [5;5]*1e3;
search_size = 5e3;
epsilon = .1e3;

x_est_true = tdoa.mlSoln(x_tdoa,zeta_true,C_roa,x_ctr,search_size,epsilon,n_tdoa);
x_est = tdoa.mlSoln(x_tdoa,zeta,C_roa,x_ctr,search_size,epsilon,n_tdoa);

fprintf('True ML Est.: (%.2f, %.2f) km, error: %.2f km\n',...
        x_est_true(1)/1e3, x_est_true(2)/1e3, norm(x_est_true-x_tgt)/1e3);
fprintf('Biased ML Est.: (%.2f, %.2f) km, error: %.2f km\n',...
        x_est(1)/1e3, x_est(2)/1e3, norm(x_est-x_tgt)/1e3);

%% ML Solvers with Bias Estimation
th_ctr = cat(1,x_ctr,zeros(n_tdoa,1),x_tdoa(:)); % expand parameter vector
search_size = cat(1,5e3*ones(2,1), 80*ones(n_tdoa-1,1), 0, zeros(numel(x_tdoa),1)); % search +/- 10 km
epsilon = cat(1,500*ones(2,1), 10*ones(n_tdoa,1), ones(numel(x_tdoa),1));
C_beta = .001*eye(numel(x_tdoa)); % position covariance error -- doesn't really matter, since we're not searching over it

[x_est_bias, alpha_est, ~] = tdoa.mlSolnUnc(x_tdoa, zeta, C_roa, C_beta,th_ctr, search_size, epsilon, n_tdoa);
fprintf('ML Est. w/Uncertainty: (%.2f, %.2f) km, error: %.2f km\n',...
        x_est_bias(1)/1e3, x_est_bias(2)/1e3, norm(x_est_bias-x_tgt)/1e3);

fprintf('True range bias: (%d, %d, %d, %d) m\n', alpha_tdoa);
fprintf('Estimated range bias: (%d, %d, %d, %d) m\n', alpha_est(1:n_tdoa));

%% Plots Solutions
fig3 = figure;
[C,h]=contourf(x_vec_km,x_vec_km,ell_true_plot,levels);
clabel(C,h);
utils.excludeFromLegend(h);
hold on;set(gca,'ydir','normal');
scatter(x_tdoa(1,:)/1e3,x_tdoa(2,:)/1e3,'^','filled','DisplayName','Sensors');
scatter(x_tgt(1)/1e3,x_tgt(2)/1e3,'o','filled','DisplayName','Target');
scatter(x_est_true(1)/1e3,x_est_true(2)/1e3,'v','filled','DisplayName','ML Est.')
legend('Location','NorthWest');
colorbar;
colormap(utils.viridis);
caxis([-100,0]);
legend('Location','NorthWest');
grid on;

utils.setPlotStyle(gca,{'tight'});

fig4 = figure;
[C,h]=contourf(x_vec_km,x_vec_km,ell_plot,levels);
clabel(C,h);
utils.excludeFromLegend(h);
hold on;set(gca,'ydir','normal');
scatter(x_tdoa(1,:)/1e3,x_tdoa(2,:)/1e3,'^','filled','DisplayName','Sensors');
scatter(x_tgt(1)/1e3,x_tgt(2)/1e3,'o','filled','DisplayName','Target');
scatter(x_est(1)/1e3,x_est(2)/1e3,'v','filled','DisplayName','ML Est.')
scatter(x_est_bias(1)/1e3,x_est_bias(2)/1e3,'s','filled','DisplayName','ML Est. w/uncertainty');
legend('Location','NorthWest');
colorbar;
colormap(utils.viridis);
caxis([-100,0]);
grid on;
legend('Location','NorthWest');

utils.setPlotStyle(gca,{'tight'});




figs = [fig1, fig2, fig3, fig4];