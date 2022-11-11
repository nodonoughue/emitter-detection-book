function figs = book2_ex3_4()
% figs=book2_ex3_4()
%
% Executes Example 3.4 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle
%
% Nicholas O'Donoughue
% 4 June 2021

% Turn off ill conditioned matrix warnings
warning('off','MATLAB:nearlySingularMatrix');

% Set up sensor and target coordinates
num_mc = 1000;
x_source_ctr = [3; 4]*1e3;
x_source = x_source_ctr + 1e2*(-1 + 2*rand(2,num_mc));

x_tdoa = [1, 3, 4, 5, 2; 0, .5, 0, .5, -1]*1e3;
n_tdoa = size(x_tdoa,2);

% Error Covariance Matrix
err_time = 1e-7; % 100 ns timing error
cov_t = err_time^2 * eye(n_tdoa);   % Time-of-Arrival
cov_r = cov_t * utils.constants.c^2;  % Range-of-Arrival

% TDOA measurement and combined covariance matrix
z_cmn = tdoa.measurement(x_tdoa, x_source);  % RDOA
z_full = tdoa.measurement(x_tdoa, x_source, 'full'); % RDOA

[test_tdoa,ref_tdoa] = utils.parseReferenceSensor(1,n_tdoa);
cov_z = utils.resampleCovMtx(cov_r, test_tdoa, ref_tdoa);
cov_z = utils.ensureInvertible(cov_z);
[test_tdoa,ref_tdoa] = utils.parseReferenceSensor('full',n_tdoa);
cov_z_full = utils.resampleCovMtx(cov_r, test_tdoa, ref_tdoa);
cov_z_full = utils.ensureInvertible(cov_z_full);

% Generate Random Noise
noise_sensor = chol(cov_r,'lower')*randn(n_tdoa, num_mc);
noise = utils.resampleNoise(noise_sensor, 1, []);
noise_full = utils.resampleNoise(noise_sensor, 'full', []);

% Noisy Measurements
zeta = z_cmn + noise;
zeta_full = z_full + noise_full;

%% ML Search Parameters
x_ctr = [2.5; 2.5]*1e3;
grid_size = [5e3; 5e3];
grid_res = 20;  % meters, grid resolution

%% GD and LS Search Parameters
x_init = [1; 1]*1e3;
epsilon = grid_res;
max_num_iterations = 200;
force_full_calc = true;
plot_progress = false;

rmse_ml = zeros(num_mc, 1);
rmse_gd = zeros(num_mc, max_num_iterations);
rmse_ls = zeros(num_mc, max_num_iterations);

rmse_ml_full = zeros(num_mc, 1);
rmse_gd_full = zeros(num_mc, max_num_iterations);
rmse_ls_full = zeros(num_mc, max_num_iterations);

for idx=1:num_mc
    fprintf('.');
    if mod(idx,100)==0
        fprintf('\n');
    end
    
    %% ML Soln
    x_ml = tdoa.mlSoln(x_tdoa, zeta(:,idx), cov_r, x_ctr, grid_size, epsilon);
    x_ml_full = tdoa.mlSoln(x_tdoa, zeta_full(:,idx), cov_r, x_ctr, grid_size, epsilon,'full');


    %% GD Soln
    [x_gd, x_gd_iters] = tdoa.gdSoln(x_tdoa, zeta(:,idx), cov_r, x_init,... 
                         [], [], epsilon, max_num_iterations, force_full_calc, plot_progress);
    [x_gd_full, x_gd_full_iters] = tdoa.gdSoln(x_tdoa, zeta_full(:,idx), cov_r, x_init,... 
                         [], [], epsilon, max_num_iterations, force_full_calc, plot_progress, 'full');

                     
    %% LS Soln
    [x_ls, x_ls_iters] = tdoa.lsSoln(x_tdoa, zeta(:,idx), cov_r, x_init,...
                         epsilon, max_num_iterations, force_full_calc, plot_progress);
    [x_ls_full, x_ls_full_iters] = tdoa.lsSoln(x_tdoa, zeta_full(:,idx), cov_r, x_init,...
                         epsilon, max_num_iterations, force_full_calc, plot_progress, 'full');

    rmse_ml(idx) = norm(x_ml-x_source(:,idx));
    rmse_gd(idx,:) = sqrt(sum(abs(x_gd_iters-x_source(:,idx)).^2,1));
    rmse_ls(idx,:) = sqrt(sum(abs(x_ls_iters-x_source(:,idx)).^2,1));

    rmse_ml_full(idx) = norm(x_ml_full-x_source(:,idx));
    rmse_gd_full(idx,:) = sqrt(sum(abs(x_gd_full_iters-x_source(:,idx)).^2,1));
    rmse_ls_full(idx,:) = sqrt(sum(abs(x_ls_full_iters-x_source(:,idx)).^2,1));
end

rmse_avg_ml = sum(rmse_ml)/num_mc;
rmse_avg_gd = sum(rmse_gd,1)/num_mc;
rmse_avg_ls = sum(rmse_ls,1)/num_mc;

rmse_avg_ml_full = sum(rmse_ml_full)/num_mc;
rmse_avg_gd_full = sum(rmse_gd_full,1)/num_mc;
rmse_avg_ls_full = sum(rmse_ls_full,1)/num_mc;

fig1=figure;
plot([1 max_num_iterations],rmse_avg_ml*[1 1],'DisplayName','ML');
hold on;
plot(1:max_num_iterations, rmse_avg_gd,'DisplayName','Gradient Descent');
plot(1:max_num_iterations, rmse_avg_ls,'DisplayName','Least Squares')
set(gca,'ColorOrderIndex',1);
plot([1 max_num_iterations],rmse_avg_ml_full*[1 1],'--','DisplayName','ML (full)');
plot(1:max_num_iterations, rmse_avg_gd_full,'--','DisplayName','Gradient Descent (full)');
plot(1:max_num_iterations, rmse_avg_ls_full,'--','DisplayName','Least Squares (full)')
set(gca,'yscale','log');
legend('Location','NorthEast')


%% Generate CRLB
% Compute the CRLB
crlb = tdoa.computeCRLB(x_tdoa, x_source_ctr, cov_t);
crlb_full = tdoa.computeCRLB(x_tdoa, x_source_ctr, cov_t, 'full');
display(crlb);
display(crlb_full);

% Compute and display the RMSE
rmse_crlb =sqrt(trace(crlb));
rmse_crlb_full =sqrt(trace(crlb_full));

fprintf('RMSE: %.2f km (%.2f km using full set)\n',rmse_crlb/1e3, rmse_crlb_full/1e3);
plot([1 max_num_iterations], rmse_crlb*[1 1],'k','DisplayName','CRLB');
plot([1 max_num_iterations], rmse_crlb_full*[1 1],'k','DisplayName','CRLB (full)');
xlabel('Iteration Number');
ylabel('RMSE [m]');
title('Monte Carlo Geolocation Results');
xlim([1 max_num_iterations]);

% Compute and display the CEP50
cep50 = utils.computeCEP50(crlb);
cep50_full = utils.computeCEP50(crlb_full);
fprintf('CEP50: %.2f km (%.2f km using full set)\n',cep50/1e3, cep50_full/1e3);

% Generate the 90% error ellipse from the CRLB
crlb_ellipse = utils.drawErrorEllipse(x_source(:,end), crlb, 101, 50);

%% Plot Result

fig2=figure;
plot(x_source(1, end), x_source(2, end), 'kx', 'DisplayName','Target');
hold on;
plot(x_tdoa(1, :), x_tdoa(2, :), 'ks', 'DisplayName','TDOA Sensor');

plot(x_ml(1), x_ml(2), 'v', 'DisplayName', 'ML Solution');
plot(x_ml_full(1), x_ml_full(2), '^', 'DisplayName', 'ML Solution (full)');
hdl=plot(x_gd_iters(1,:), x_gd_iters(2,:), '-.');
utils.excludeFromLegend(hdl);
plot(x_gd(1),x_gd(2),'-.+','DisplayName','GD Solution','Color',hdl.Color);
hdl=plot(x_gd_full_iters(1,:), x_gd_full_iters(2,:), '-.');
utils.excludeFromLegend(hdl);
plot(x_gd_full(1),x_gd_full(2),'-.+','DisplayName','GD Solution (full)','Color',hdl.Color);
hdl=plot(x_ls_iters(1,:), x_ls_iters(2,:), '-');
utils.excludeFromLegend(hdl);
plot(x_ls(1), x_ls(2), '-*','DisplayName','LS Solution','Color',hdl.Color);
hdl=plot(x_ls_full_iters(1,:), x_ls_full_iters(2,:), '-');
utils.excludeFromLegend(hdl);
plot(x_ls_full(1), x_ls_full(2), '-*','DisplayName','LS Solution (full)','Color',hdl.Color);

plot(crlb_ellipse(1,:), crlb_ellipse(2,:), '--','DisplayName','CRLB');

grid on;
ylim([-1 6]*1e3);
xlim([0.5 5.5]*1e3);
caxis([-20 0]);
set(gca,'ydir','normal');
legend('Location','NorthEast');
utils.setPlotStyle(gca,{'widescreen'});

%% Turn on ill conditioned matrix warnings
warning('on','MATLAB:nearlySingularMatrix');

%% Parse Outputs
figs = [fig1, fig2];