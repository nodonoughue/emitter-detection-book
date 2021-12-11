function figs = book2_ex4_1()
% figs=book2_ex4_1()
%
% Executes Example 4.1 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle
%
% Nicholas O'Donoughue
% 24 June 2021

% Set up sensors
x_sensor_enu = [-5, 0, 5;
                 0, 0, 0;
                 0, 0, 0]*1e3;
kph2mps = 1e3/60^2;
v_sensor_enu = [80, 80, 80;
                0, 0, 0;
                0, 0, 0]*kph2mps;  % Convert to SI units (m/s)
ref_lat = 5; % deg Lat (N)
ref_lon = -15; % deg Lon (W)
ref_alt = 10e3; % alt (m)

err_aoa_deg = 2;
err_tdoa_s = 1e-6;
err_fdoa_Hz = 10;

x_source_lla = [4.5, -(14 + 40/60), 0]'; % deg Lat (N), deg Lon (E), alt (m)
v_source_enu = [0, 0, 0];
f_source_Hz = 1e9;

%% Convert positions and velocity to ECEF
% [x, y, z] = utils.enu2ecef(x_sensor_enu(1,:), x_sensor_enu(2,:), x_sensor_enu(3,:), ...
%                            ref_lat, ref_lon, ref_alt);
% x_sensor_ecef = [x(:), y(:), z(:)]';
% 
% [vx, vy, vz] = utils.enu2ecef_vel(v_sensor_enu(1,:), v_sensor_enu(2,:), v_sensor_enu(3,:),...
%                                   ref_lat, ref_lon);
% v_sensor_ecef = [vx(:), vy(:), vz(:)]';
% 
% [x, y, z] = utils.lla2ecef(x_source_lla(1), x_source_lla(2), ...
%                            x_source_lla(3), 'deg');
% x_source_ecef = [x(:), y(:), z(:)]';
% v_source_ecef = zeros(size(x_source_ecef)); % stationary emitter

[e_source, n_source, u_source] = utils.lla2enu(x_source_lla(1), x_source_lla(2), x_source_lla(3),...
                                               ref_lat, ref_lon, ref_alt);
x_source_enu = [e_source; n_source; u_source];

%% Sensor Selection
n_tdoa = size(x_sensor_enu,2);
n_aoa = n_tdoa;
n_fdoa = n_aoa;
ref_tdoa = n_tdoa;
ref_fdoa = n_fdoa;

% Manually do a reference/text index set
% 1:n_aoa are AOA sensors
% n_aoa + (1:n_tdoa) are TDOA sensors
% n_aoa + n_tdoa + (1:n_fdoa) are FDOA sensors
[tdoa_ref_vec, tdoa_test_vec] = utils.parseReferenceSensor(ref_tdoa, n_tdoa);
[fdoa_ref_vec, fdoa_test_vec] = utils.parseReferenceSensor(ref_fdoa, n_fdoa);
ref_idx =  [1:2*n_aoa,      2*n_aoa + tdoa_ref_vec,  2*n_aoa + n_tdoa + fdoa_ref_vec];
test_idx = [nan(1,2*n_aoa), 2*n_aoa + tdoa_test_vec, 2*n_aoa + n_tdoa + fdoa_test_vec];

%% Error Covariance Matrix
cov_psi = (err_aoa_deg*pi/180)^2*eye(2*n_aoa); % rad^2
cov_r = (err_tdoa_s*utils.constants.c)^2*eye(n_tdoa); % m^2
cov_rr = (err_fdoa_Hz*utils.constants.c/f_source_Hz)^2*eye(n_fdoa); % m^2/s^2

cov_x = blkdiag(cov_psi, cov_r, cov_rr);

%% Generate Noise
num_mc = 10000;
noise_white = randn(2*n_aoa+n_tdoa+n_fdoa, num_mc); % one row per sensor,
                                                    % one column per trial
% Generate sensor level (n_aoa + n_tdoa + n_fdoa) noise with proper errors
noise_sensor = chol(cov_x,'lower')*noise_white;
% Resample to account for reference sensors used in TDOA and FDOA
noise_msmt = utils.resampleNoise(noise_sensor,test_idx, ref_idx);

%% Generate Data
z = hybrid.measurement(x_sensor_enu, x_sensor_enu, x_sensor_enu, ...
                       v_sensor_enu - v_source_enu, x_source_enu, ref_tdoa, ref_fdoa);
zeta = z + noise_msmt;

%% GD and LS Search Parameters
x_init_enu = [1; 1; 0]*1e3;
% [xx, yy, zz] = utils.enu2ecef(x_init_enu(1), x_init_enu(2), x_init_enu(3), ref_lat, ref_lon, ref_alt);
% x_init_ecef = [xx;yy;zz]; % reassemble ecef coordinates into a single variable

epsilon = 100;
max_num_iterations = 100;
force_full_calc = true;
plot_progress = false;

error = zeros(num_mc, max_num_iterations);

for idx=1:num_mc
    if mod(idx,10)==0
        fprintf('.');
    end
    if mod(idx,1000)==0
        fprintf(' (%d/%d)\n', idx, num_mc);
    end
    
    %% TDOA, AOA, and FDOA Error
    
    %% LS Soln
    [x_ls, x_ls_iters] = hybrid.lsSoln(x_sensor_enu, x_sensor_enu, x_sensor_enu, ...
                                       v_sensor_enu, zeta(:,idx), cov_x, x_init_enu,...
                                       epsilon, max_num_iterations, force_full_calc, plot_progress,...
                                       ref_tdoa, ref_fdoa);

    error(idx,:) = sqrt(sum(abs(x_ls_iters-x_source_enu).^2,1));

end

% Remove outliers
error = rmoutliers(error);

%% Plot RMSE and CRLB
rmse_ls = sqrt(sum(error.^2,1)/num_mc);

fig2=figure;
plot(1:max_num_iterations, rmse_ls,'DisplayName','Least Squares')
set(gca,'yscale','log');
legend('Location','NorthEast')

% Compute the CRLB
crlb = hybrid.computeCRLB(x_sensor_enu, x_sensor_enu, x_sensor_enu, v_sensor_enu, ...
                          x_source_enu, cov_x);

% Compute and display the RMSE
rmse_crlb =sqrt(trace(crlb));

hold on;
plot([1 max_num_iterations], rmse_crlb*[1 1],'k','DisplayName','CRLB');
xlabel('Iteration Number');
ylabel('RMSE [m]');
title('Monte Carlo Geolocation Results');

%% Compute and display the CEP50
cep50 = utils.computeCEP50(crlb);
fprintf('CEP50: %.2f km \n',cep50/1e3);

% Generate the 90% error ellipse from the CRLB
crlb_ellipse = utils.drawErrorEllipse(x_source_enu(1:2), crlb, 101, 50);

%% Plot Result

fig1=figure;
plot(x_source_enu(1, end), x_source_enu(2, end), 'kx', 'DisplayName','Target');
hold on;
plot(x_sensor_enu(1, :), x_sensor_enu(2, :), 'ks', 'DisplayName','Sensor');

hdl=plot(x_ls_iters(1,:), x_ls_iters(2,:), '-');
utils.excludeFromLegend(hdl);
plot(x_ls(1), x_ls(2), '-*','DisplayName','LS Solution','Color',hdl.Color);

plot(crlb_ellipse(1,:), crlb_ellipse(2,:), '--','DisplayName','CRLB');

grid on;
set(gca,'ydir','normal');
legend('Location','NorthEast');

%% Parse Outputs
figs = [fig1, fig2];