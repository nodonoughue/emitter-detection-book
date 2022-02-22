function figs = book2_ex6_5()
% figs = book2_ex6_5()
%
% Executes Example 6.5 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   figs         array of figure handle
%
% Nicholas O'Donoughue
% 18 February 2022

%% Set up sensors
x_tdoa = [-1, 0, 1;
          0, 1, 0]*1e3;
x_fdoa = x_tdoa;
v_fdoa = [0, 0, 0;
          500, 500, 500];

n_dim = size(x_fdoa,1);
n_fdoa = size(x_fdoa,2);
n_tdoa = size(x_tdoa,2);

C_vel = 100^2*eye(n_dim*n_fdoa);

%% Generate Random Velocity Errors
U_vel = chol(C_vel,'lower');
vel_err = reshape(U_vel * randn(n_dim*n_fdoa,1),n_dim,n_fdoa);

v_fdoa_actual = v_fdoa + vel_err;

%% Generate Measurements
x_tgt = [-3;4] * 1e3;
x_cal = [-2,-1,0,1,2;-5,-5,-5,-5,-5] * 1e3;
n_cal = size(x_cal,2);

z = hybrid.measurement([],x_tdoa,x_fdoa,v_fdoa_actual,x_tgt, n_tdoa, n_fdoa);
z_cal = hybrid.measurement([],x_tdoa,x_fdoa,v_fdoa_actual,x_cal, n_tdoa,n_fdoa); % free of pos unc and bias

% Build sensor-level covariance matrix
err_toa = 100e-9;
err_foa = 1;
f0 = 10e9;
lam = utils.constants.c/f0;
C_toa = err_toa^2*eye(n_tdoa);
C_roa = utils.constants.c^2*C_toa; % range error covariance
C_foa = err_foa^2*eye(n_fdoa);
C_rroa = lam^2*C_foa; % range rate error covariance
C_tf = blkdiag(C_roa, C_rroa); % assume TDOA/FDOA measurements are independent

% Generate noise; we'll do it separately for TDOA and FDOA
C_rdoa = utils.resampleCovMtx(C_roa,n_tdoa);
L_rdoa= chol(C_rdoa,'lower');
noise_rdoa = L_rdoa*randn(size(L_rdoa,2),n_cal+1);

C_rrdoa = utils.resampleCovMtx(C_rroa,n_tdoa);
L_rrdoa= chol(C_rrdoa,'lower');
noise_rrdoa = L_rrdoa*randn(size(L_rrdoa,2),n_cal+1);

noise = cat(1,noise_rdoa,noise_rrdoa);
zeta = z + noise(:,1);
zeta_cal = z_cal + noise(:,2:end);


%% Estimate Position
x_init = [0; 5]*1e3;

[x_est,x_est_full] = hybrid.gdSoln([],x_tdoa,x_fdoa,v_fdoa,zeta,C_tf,x_init,[],[],[],[],false,false,n_tdoa,n_fdoa);
[x_est_cal,x_est_cal_full,~,~] = hybrid.gdSolnCal([],x_tdoa,x_fdoa,v_fdoa,zeta,x_cal,zeta_cal,C_tf,x_init,[],[],[],[],false,false,n_tdoa,n_fdoa);


%% Plot Scenario
fig1 = figure;
hold on;
scatter(x_tdoa(1,:)/1e3,x_tdoa(2,:)/1e3,'^','filled','DisplayName','Sensors');
scatter(x_tgt(1)/1e3,x_tgt(2)/1e3,'o','filled','DisplayName','Target');
scatter(x_cal(1,:)/1e3,x_cal(2,:)/1e3,'v','filled','DisplayName','Calibration Emitter');
scatter(x_init(1)/1e3,x_init(2)/1e3,'+','DisplayName','Initial Guess');
hdl = plot(x_est_full(1,:)/1e3,x_est_full(2,:)/1e3,'--');
utils.excludeFromLegend(hdl);
scatter(x_est(1)/1e3,x_est(2)/1e3,'s','filled','MarkerFaceColor',hdl.Color,'DisplayName','Solution (w/o cal)');
hdl = plot(x_est_cal_full(1,:)/1e3,x_est_cal_full(2,:)/1e3,'--');
utils.excludeFromLegend(hdl);
scatter(x_est_cal(1)/1e3,x_est_cal(2)/1e3,'s','filled','MarkerFaceColor',hdl.Color,'DisplayName','Solution (w/cal)');
grid on
legend('Location','NorthEast');

utils.setPlotStyle(gca,{'tight'});


%%  Bonus: FDOA-only Cal
z_cal2 = fdoa.measurement(x_fdoa,v_fdoa_actual,x_cal, n_fdoa); % free of pos unc and bias
zeta_cal2 = z_cal2 + noise_rrdoa(:,2:end);
n_cal_msmt = numel(zeta_cal2);


% Manually estimate velocity error term
y_b = @(beta) reshape(zeta_cal2 - fdoa.measurement(reshape(beta(1:n_dim*n_fdoa),n_dim,n_fdoa), ...
                                                   reshape(beta(n_dim*n_fdoa+(1:n_dim*n_fdoa)),n_dim,n_fdoa), ... % v_fdoa
                                                   x_cal, n_fdoa), ...
                                                   n_cal_msmt,1);

J_b = @(beta) reshape(fdoa.grad_b(reshape(beta(1:n_dim*n_fdoa),n_dim,n_fdoa), ...
                                    reshape(beta(n_dim*n_fdoa+(1:n_dim*n_fdoa)),n_dim,n_fdoa), ... % v_fdoa
                                    x_cal, n_fdoa), ...
                                    n_dim*2*n_fdoa,n_cal_msmt);

% Build the initial vector
beta_init = [x_fdoa(:); v_fdoa(:)];
C_cal = kron(eye(n_cal),C_rrdoa);
beta_est = utils.gdSoln(y_b,J_b,C_cal,beta_init,[],[],[],[],false,false);
v_fdoa_est = reshape(beta_est(n_dim*n_fdoa + (1:n_dim*n_fdoa)),n_dim,n_fdoa);

[x_est_cal2,x_est_cal2_full] = hybrid.gdSoln([],x_tdoa,x_fdoa,v_fdoa_est,zeta,C_tf,x_init,[],[],[],[],false,false,n_tdoa,n_fdoa);

fig2 = figure;
hold on;
scatter(x_tdoa(1,:)/1e3,x_tdoa(2,:)/1e3,'^','filled','DisplayName','Sensors');
scatter(x_tgt(1)/1e3,x_tgt(2)/1e3,'o','filled','DisplayName','Target');
scatter(x_cal(1,:)/1e3,x_cal(2,:)/1e3,'v','filled','DisplayName','Calibration Emitter');
scatter(x_init(1)/1e3,x_init(2)/1e3,'+','DisplayName','Initial Guess');
hdl = plot(x_est_full(1,:)/1e3,x_est_full(2,:)/1e3,'--');
utils.excludeFromLegend(hdl);
scatter(x_est(1)/1e3,x_est(2)/1e3,'s','filled','MarkerFaceColor',hdl.Color,'DisplayName','Solution (w/o cal)');
hdl = plot(x_est_cal_full(1,:)/1e3,x_est_cal_full(2,:)/1e3,'--');
utils.excludeFromLegend(hdl);
scatter(x_est_cal(1)/1e3,x_est_cal(2)/1e3,'s','filled','MarkerFaceColor',hdl.Color,'DisplayName','Solution (w/cal)');
hdl = plot(x_est_cal2_full(1,:)/1e3,x_est_cal2_full(2,:)/1e3,'--');
utils.excludeFromLegend(hdl);
scatter(x_est_cal2(1)/1e3,x_est_cal2(2)/1e3,'s','filled','MarkerFaceColor',hdl.Color,'DisplayName','Solution (w/FDOA cal)');
grid on
legend('Location','NorthEast');
xlim([-4,4]);
ylim([-6,10]);
utils.setPlotStyle(gca,{'tight'});

figs = [fig1,fig2];