function figs = book2_ex5_5()
% fig=book2_ex5_5()
%
% Executes Example 5.5 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle
%
% Nicholas O'Donoughue
% 24 November 2021

%% Set up scene
% Positions
baseline = 10e3;
n_tdoa = 4;
tdoa_ang = pi/6 + 2*pi/3 *(0:n_tdoa-2); % for the outer stations in the Y configuration
x_tdoa = baseline * [0, cos(tdoa_ang); ...
                     0, sin(tdoa_ang);
                     0, 0, 0, 0]; % ENU

time_err = 3e-7;
C_roa = (utils.constants.c*time_err)^2 * eye(n_tdoa);

% Define target coordinates
tgt_rng = 100e3;
tgt_alt = 40e3*unitsratio('m','ft');
x_tgt_g = [tgt_rng; 0; tgt_alt]; % spherical Earth coordinates (East, North, Alt)
[e, n, u] = utils.correctENU(x_tgt_g(1), x_tgt_g(2), x_tgt_g(3));
x_tgt = [e; n; u]; % ENU

%% External Prior
x_prior_g = [95; 10; 10]*1e3;
[e, n, u] = utils.correctENU(x_prior_g(1), x_prior_g(2), x_prior_g(3));
x_prior = [e; n; u];
C_prior = [5, 1, 0; 1, 50, 0; 0, 0, 10]*1e6;

% Note that MVNPDF wants to work on rows, not columns; so let's transpose
% x and x_prior
prior = @(x) utils.mvnpdf(x', x_prior', C_prior);

%% Measurement
z = tdoa.measurement(x_tdoa, x_tgt);
C_rdoa = utils.resampleCovMtx(C_roa,[]);
L = chol(C_rdoa,'lower');
n = L*randn(size(z));

zeta = z+n;

%% Solution
x_ctr = x_tgt;
grid_size = [50e3, 50e3, 0];
epsilon = 250;

[x_ml, A, x_grid] = tdoa.mlSoln(x_tdoa, zeta, C_roa, x_ctr, grid_size,epsilon);
[x_ml_p, A_p, ~] = tdoa.mlSolnPrior(x_tdoa, zeta, C_roa, prior, x_ctr, grid_size, epsilon);

fprintf('Solution w/o prior: %.2f km, %.2f km, %.2f km\n',x_ml/1e3);
fprintf('\tError: %.2f km\n', norm(x_ml-x_tgt)/1e3);
fprintf('Solution w/prior:   %.2f km, %.2f km, %.2f km\n',x_ml_p/1e3);
fprintf('\tError: %.2f km\n', norm(x_ml_p-x_tgt)/1e3);

%% Plot
x_vec = x_grid{1};
y_vec = x_grid{2};

fig1=figure;
imagesc(x_vec, y_vec, reshape(A,numel(x_vec),numel(y_vec))');
set(gca,'ydir','normal');
hold on;
scatter(x_tdoa(1,:),x_tdoa(2,:),'o','filled','DisplayName','Sensors')
scatter(x_tgt(1), x_tgt(2),'^','filled','DisplayName','Target');
scatter(x_ml(1), x_ml(2), 's','filled','DisplayName','Estimate');

grid on;
legend('Location','NorthWest');
colorbar;
caxis([-50 0]);
colormap(utils.viridis);

xlim([min(x_tdoa(1,:))-10e3, max(x_vec)]);
xlabel('x [m]');
ylabel('y [m]');
title('Likelihood and Estimate without Prior');

fig2=figure;
imagesc(x_vec, y_vec, reshape(A_p, numel(x_vec), numel(y_vec))');
set(gca,'ydir','normal');
hold on;
scatter(x_tdoa(1,:),x_tdoa(2,:),'o','filled','DisplayName','Sensors')
scatter(x_tgt(1), x_tgt(2),'^','filled','DisplayName','Target');
scatter(x_ml(1), x_ml(2), 's','filled','DisplayName','Estimate (w/o prior)');
scatter(x_prior(1), x_prior(2), 'v','filled','DisplayName','Prior');
ell = utils.drawErrorEllipse(x_prior(1:2),C_prior(1:2,1:2),101,90);
plot(ell(1,:),ell(2,:),'-.','Color','w','DisplayName','Prior Confidence (90%)')
scatter(x_ml_p(1), x_ml_p(2), 'd','filled','DisplayName','Estimate (w/prior)');

grid on;
legend('Location','NorthWest');
colorbar;
caxis([-50 0]);
colormap(utils.viridis);

xlim([min(x_tdoa(1,:))-10e3, max(x_vec)]);
xlabel('x [m]');
ylabel('y [m]');
title('Likelihood and Estimate with Prior');

%% Collect Figure Handles for Export
figs = [fig1, fig2];