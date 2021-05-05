function figs = book2_ex2_1()
% figs=book2_ex2_1()
%
% Executes Example 2.1 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   figs         array of figures
%
% Nicholas O'Donoughue
% 22 April 2021

% Set up sensor and target coordinates
x_source = [3; 3]*1e3;

x_aoa = [4; 0]*1e3;
x_tdoa = [1, 3; 0, .5]*1e3;
x_fdoa = [0, 0; 1, 2]*1e3;
v_fdoa = [1, 1; -1, -1]*sqrt(.5)*300; % 300 m/s, at -45 deg heading

% Draw the geometry
fig1=figure;

plot(x_source(1), x_source(2), 'kx', 'DisplayName','Target');
hold on;
h_aoa=plot(x_aoa(1), x_aoa(2), 'o', 'DisplayName','AOA Sensor');
h_tdoa=plot(x_tdoa(1, :), x_tdoa(2, :), 's', 'DisplayName','TDOA Sensor');
h_fdoa=plot(x_fdoa(1, :), x_fdoa(2, :), '^', 'DisplayName','FDOA Sensor');
h=utils.drawArrow(x_fdoa(1,1)+[0 v_fdoa(1,1)],x_fdoa(2,1)+[0 v_fdoa(2,1)]);
h.Color = h_fdoa.Color; % Make sure the arrow color matches the marker
h=utils.drawArrow(x_fdoa(1,2)+[0 v_fdoa(1,2)],x_fdoa(2,2)+[0 v_fdoa(2,2)]);
h.Color = h_fdoa.Color; % Make sure the arrow color matches the marker

% AOA/TDOA/FDOA Solutions
psi = triang.measurement(x_aoa, x_source, false);
aoa_soln = triang.drawLob(x_aoa,psi,x_source,1.5);
plot(aoa_soln(1,:), aoa_soln(2,:), 'Color', h_aoa.Color, 'DisplayName', 'AOA Solution');

r = tdoa.measurement(x_tdoa, x_source);
tdoa_soln = tdoa.drawIsochrone(x_tdoa(:,2),x_tdoa(:,1),r,1000,5e3);
plot(tdoa_soln(1,:), tdoa_soln(2,:), 'Color', h_tdoa.Color, 'DisplayName', 'TDOA Solution');

r_dot = fdoa.measurement(x_fdoa, v_fdoa, x_source);
soln_fdoa = fdoa.drawIsodop(x_fdoa(:,2),v_fdoa(:,2),x_fdoa(:,1),v_fdoa(:,1),r_dot,2000,15e3);
plot(soln_fdoa(1,:), soln_fdoa(2,:), 'Color', h_fdoa.Color, 'DisplayName', 'FDOA Solution');

ylim([0 4]*1e3);
xlim([-0.5 5.5]*1e3);

legend('location','NorthEast')
utils.setPlotStyle(gca,{'widescreen'});

%% Compute Variances and Print
err_aoa = 3; % deg
cov_psi = (err_aoa*pi/180)^2; % rad^2
fprintf('AOA Measurement: %.2f deg\n',psi*180/pi);
fprintf('AOA Covariance: %.2g rad^2\n',cov_psi);

err_time = 1e-7; % 100 ns timing error
err_r = err_time * utils.constants.c;
cov_r = 2 * (err_r)^2; % m^2, double for the combination of test/ref msmts
fprintf('TDOA Measurement: %.2f m\n',r);
fprintf('TDOA Covariance: %.2g m^2\n',cov_r);

freq_err = 10; % Hz
f0 = 1e9; % Hz
rr_err = freq_err * utils.constants.c/f0; % (m/s)
cov_rr = 2 * rr_err^2; % (m/s)^2
fprintf('FDOA Measurement: %.2f m/s\n',r_dot);
fprintf('FDOA Covariance: %.2g (m/s)^2\n',cov_rr);

z = hybrid.measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source);
cov_z = diag([cov_psi, cov_r, cov_rr]);

%% Set Up Search Grid
x_grid = (-.5:.02:5.5)*1e3;
y_grid = (0:.02:4)*1e3;
[xx,yy] = meshgrid(x_grid, y_grid);
x_test_pos = [xx(:), yy(:)]'; % 2 x n_test_pos

%% Plot AOA Likelihood
ell_aoa = reshape(triang.loglikelihood(x_aoa, psi, cov_psi, x_test_pos), size(xx));

fig2=figure;
imagesc(x_grid,y_grid,ell_aoa);
hold on;
colorbar;
cmap = utils.viridis;
colormap(flipud(cmap));

plot(x_source(1), x_source(2), 'wx', 'DisplayName','Target');
plot(x_aoa(1), x_aoa(2), 'ko', 'DisplayName','AOA Sensor');

ylim([0 4]*1e3);
xlim([-0.5 5.5]*1e3);
caxis([-20 0]);
set(gca,'ydir','normal');
legend('Location','NorthEast');

utils.setPlotStyle(gca,{'widescreen'});

%% TDOA
ell_tdoa = reshape(tdoa.loglikelihood(x_tdoa, r, cov_r, x_test_pos), size(xx));

fig3=figure;
imagesc(x_grid,y_grid,ell_tdoa);
hold on;
colorbar;
cmap = utils.viridis;
colormap(flipud(cmap));

plot(x_source(1), x_source(2), 'wx', 'DisplayName','Target');
plot(x_tdoa(1, :), x_tdoa(2, :), 'ks', 'DisplayName','TDOA Sensor');

ylim([0 4]*1e3);
xlim([-0.5 5.5]*1e3);
caxis([-20 0]);
set(gca,'ydir','normal');
legend('Location','NorthEast');

utils.setPlotStyle(gca,{'widescreen'});

%% FDOA
ell_fdoa = reshape(fdoa.loglikelihood(x_fdoa, v_fdoa, r_dot, cov_rr, x_test_pos), size(xx));

fig4=figure;
imagesc(x_grid,y_grid,ell_fdoa);
hold on;
colorbar;
cmap = utils.viridis;
colormap(flipud(cmap));

plot(x_source(1), x_source(2), 'wx', 'DisplayName','Target');
plot(x_fdoa(1, :), x_fdoa(2, :), 'k^', 'DisplayName','FDOA Sensor');
utils.drawArrow(x_fdoa(1,1)+[0 v_fdoa(1,1)],x_fdoa(2,1)+[0 v_fdoa(2,1)]);
utils.drawArrow(x_fdoa(1,2)+[0 v_fdoa(1,2)],x_fdoa(2,2)+[0 v_fdoa(2,2)]);

ylim([0 4]*1e3);
xlim([-0.5 5.5]*1e3);
caxis([-20 0]);
set(gca,'ydir','normal');
legend('Location','NorthEast');

utils.setPlotStyle(gca,{'widescreen'});

%% Hybrid
ell_hybrid = reshape(hybrid.loglikelihood(x_aoa, x_tdoa, x_fdoa, v_fdoa, z, cov_z, x_test_pos), size(xx));

fig5=figure;
imagesc(x_grid,y_grid,ell_hybrid);
hold on;
colorbar;
cmap = utils.viridis;
colormap(flipud(cmap));

plot(x_source(1), x_source(2), 'wx', 'DisplayName','Target');
plot(x_aoa(1), x_aoa(2), 'ko', 'DisplayName','AOA Sensor');
plot(x_tdoa(1, :), x_tdoa(2, :), 'ks', 'DisplayName','TDOA Sensor');
plot(x_fdoa(1, :), x_fdoa(2, :), 'k^', 'DisplayName','FDOA Sensor');
utils.drawArrow(x_fdoa(1,1)+[0 v_fdoa(1,1)],x_fdoa(2,1)+[0 v_fdoa(2,1)]);
utils.drawArrow(x_fdoa(1,2)+[0 v_fdoa(1,2)],x_fdoa(2,2)+[0 v_fdoa(2,2)]);

ylim([0 4]*1e3);
xlim([-0.5 5.5]*1e3);
caxis([-20 0]);
set(gca,'ydir','normal');
legend('Location','NorthEast');

utils.setPlotStyle(gca,{'widescreen'});

%% Package figure handles
figs = [fig1, fig2, fig3, fig4, fig5];