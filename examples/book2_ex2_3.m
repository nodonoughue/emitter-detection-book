function fig = book2_ex2_3()
% fig=book2_ex2_3()
%
% Executes Example 2.3 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle
%
% Nicholas O'Donoughue
% 22 April 2021

% Set up sensor and target coordinates
x_source = [3; 3]*1e3;

x_aoa = [4; 0]*1e3;
x_tdoa = [1, 3; 0, .5]*1e3;
x_fdoa = [0, 0; 1, 2]*1e3;
v_fdoa = [1, 1; -1, -1]*sqrt(.5)*300; % 300 m/s, at -45 deg heading

% Error Covariance Matrix
err_aoa = 3; % deg
cov_psi = (err_aoa*pi/180)^2; % rad^2

err_time = 1e-7; % 100 ns timing error
err_r = err_time * utils.constants.c;
cov_r = (err_r)^2*eye(size(x_tdoa,2)); % m^2, sensor measurements
cov_r_out = utils.resampleCovMtx(cov_r,1); % m^2, sensor pairs

freq_err = 10; % Hz
f0 = 1e9; % Hz
rr_err = freq_err * utils.constants.c/f0; % (m/s)
cov_rr = rr_err^2 * eye(size(x_fdoa,2)); % (m/s)^2, sensor measurements
cov_rr_out = utils.resampleCovMtx(cov_rr,1); % (m/s)^2, sensor pairs

% Hybrid measurement and combined covariance matrix
z = hybrid.measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source);
cov_z = blkdiag(cov_psi, cov_r, cov_rr);
cov_z_out = blkdiag(cov_psi, cov_r_out, cov_rr_out);

% Generate Random Noise
L = chol(cov_z_out,'lower'); % Cholesky decomposition of the covariance matrix
noise = L*randn(size(L,2), 1);

% Noisy Measurements
zeta = z + noise;

%% ML Search Parameters
x_ctr = [2.5; 2.5]*1e3;
grid_size = [5e3; 5e3];
grid_res = 25;  % meters, grid resolution

%% GD and LS Search Parameters
x_init = [1; 1]*1e3;
epsilon = grid_res;
max_num_iterations = 50;
force_full_calc = true;
plot_progress = false;
    
%% ML Soln
x_ml = hybrid.mlSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, zeta, cov_z, x_ctr,...
    grid_size, epsilon);

%% GD Soln
[x_gd, x_gd_full] = hybrid.gdSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, zeta, cov_z, x_init,... 
                     [], [], epsilon, max_num_iterations, force_full_calc, plot_progress);

%% LS Soln
[x_ls, x_ls_full] = hybrid.lsSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, zeta, cov_z, x_init,...
                     epsilon, max_num_iterations, force_full_calc, plot_progress);

%% Generate CRLB
% Compute the CRLB
crlb = hybrid.computeCRLB(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source, cov_z);
display(crlb);

% Compute and display the RMSE
rmse_crlb =sqrt(trace(crlb));
fprintf('RMSE: %.2f km\n',rmse_crlb/1e3);

% Compute and display the CEP50
cep50 = utils.computeCEP50(crlb);
fprintf('CEP50: %.2f km\n',cep50/1e3);

% Generate the 90% error ellipse from the CRLB
crlb_ellipse = utils.drawErrorEllipse(x_source, crlb, 101, 90);

%% Plot Result

fig=figure;
plot(x_source(1), x_source(2), 'kx', 'DisplayName','Target');
hold on;
plot(x_aoa(1), x_aoa(2), 'ko', 'DisplayName','AOA Sensor');
plot(x_tdoa(1, :), x_tdoa(2, :), 'ks', 'DisplayName','TDOA Sensor');
plot(x_fdoa(1, :), x_fdoa(2, :), 'k^', 'DisplayName','FDOA Sensor');
utils.drawArrow(x_fdoa(1,1)+[0 v_fdoa(1,1)],x_fdoa(2,1)+[0 v_fdoa(2,1)]);
utils.drawArrow(x_fdoa(1,2)+[0 v_fdoa(1,2)],x_fdoa(2,2)+[0 v_fdoa(2,2)]);

plot(x_ml(1), x_ml(2), 'v', 'DisplayName', 'ML Solution');
hdl=plot(x_gd_full(1,:), x_gd_full(2,:), '-.');
utils.excludeFromLegend(hdl);
plot(x_gd(1),x_gd(2),'-.+','DisplayName','GD Solution','Color',hdl.Color);
hdl=plot(x_ls_full(1,:), x_ls_full(2,:), '-');
utils.excludeFromLegend(hdl);
plot(x_ls(1), x_ls(2), '-*','DisplayName','LS Solution','Color',hdl.Color);

plot(crlb_ellipse(1,:), crlb_ellipse(2,:), '--','DisplayName','CRLB');

grid on;
ylim([0 4]*1e3);
xlim([-0.5 5.5]*1e3);
caxis([-20 0]);
set(gca,'ydir','normal');
legend('Location','NorthEast');
utils.setPlotStyle(gca,{'widescreen'});
