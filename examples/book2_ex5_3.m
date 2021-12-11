function figs = book2_ex5_3()
% fig=book2_ex5_3()
%
% Executes Example 5.3 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle
%
% Nicholas O'Donoughue
% 17 November 2021

%% Set up scene
% Positions
ref_lla = [20, -150, 0];
x_aoa = [0, 0, 0]'; % AOA sensor at reference LLA; origin in ENU
x_tdoa = [20e3, 25e3; 0, 0; 0, 0];
n_aoa = size(x_aoa,2); n_tdoa = size(x_tdoa,2);

bng_tgt = 30; % deg E of N
tgt_rng = 50e3;
tgt_alt = 10e3;
x_tgt = [tgt_rng * sind(bng_tgt);
         tgt_rng * cosd(bng_tgt);
         tgt_alt];

% Error
err_aoa = 3*pi/180; % 1 deg
err_toa = 1e-6;
C_aoa = err_aoa^2*eye(2); % the aoa sensor will be 2D, unless we specify az-only
C_roa = (utils.constants.c*err_toa)^2*eye(n_tdoa);
C_full = blkdiag(C_aoa, C_roa);

%% CRLB Computation

C_raw = hybrid.computeCRLB(x_aoa, x_tdoa, [], [], x_tgt, C_full);

[~, a_grad] = utils.constraints.fixedAlt(tgt_alt, 'flat');
C_fix = hybrid.computeCRLBfixed(x_aoa, x_tdoa, [], [], x_tgt, C_full, a_grad);

fprintf('CRLB (unconstrained):\n');
disp(C_raw);
fprintf('CRLB (constrained):\n');
disp(C_fix);

%% Plot for x/y Grid
xy_vec = linspace(-10e3, 10e3, 201);
[X,Y] = meshgrid(xy_vec);
x_grid = [X(:), Y(:), zeros(numel(X),1)]' + x_tgt;

% Define target position and initial estimate
C_raw_grid = hybrid.computeCRLB(x_aoa, x_tdoa, [], [], x_grid, C_full);
C_fix_grid = hybrid.computeCRLBfixed(x_aoa, x_tdoa, [], [], x_grid, C_full, a_grad);

% Compute RMSE
rmse_raw = reshape(sqrt(C_raw_grid(1,1,:)+C_raw_grid(2,2,:)+C_raw_grid(3,3,:)),size(X))/1e3;
rmse_fix = reshape(sqrt(C_fix_grid(1,1,:)+C_fix_grid(2,2,:)+C_fix_grid(3,3,:)),size(X))/1e3;

% Plot
fig1=figure;
subplot(1,2,1);
[C,hdl] = contourf(xy_vec/1e3+x_tgt(1)/1e3, xy_vec/1e3+x_tgt(2)/1e3, rmse_raw);
clabel(C,hdl);
utils.excludeFromLegend(hdl);
hold on;
% scatter(x_tdoa(1,:), x_tdoa(2,:), 'o', 'filled', 'DisplayName','Sensors');
scatter(x_tgt(1)/1e3, x_tgt(2)/1e3, 'k^','filled','DisplayName','Target');
title('Unconstrained RMSE [km]');
xlabel('E [km]');
ylabel('N [km]');
colorbar;
caxis([0 20]);
colormap(flipud(utils.viridis));
legend('Location','NorthWest');

subplot(1,2,2);
[C,hdl]=contourf(xy_vec/1e3+x_tgt(1)/1e3, xy_vec/1e3+x_tgt(2)/1e3, rmse_fix);
clabel(C,hdl);
utils.excludeFromLegend(hdl);
hold on;
% scatter(x_tdoa(1,:), x_tdoa(2,:), 'o', 'filled', 'DisplayName','Sensors');
scatter(x_tgt(1)/1e3, x_tgt(2)/1e3, 'k^','filled','DisplayName','Target');
title('Constrained RMSE [km]');
xlabel('E [km]');
ylabel('N [km]');
colorbar;
caxis([0 20]);
colormap(flipud(utils.viridis));
legend('Location','NorthWest');

%% Collect Figure Handles for Export
figs = [fig1];