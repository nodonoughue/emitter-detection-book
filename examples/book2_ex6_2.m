function fig = book2_ex6_2()
% fig=book2_ex6_2()
%
% Executes Example 6.2 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle
%
% Nicholas O'Donoughue
% 28 November 2021

%% Set up sensors
x_tdoa = [0, 2, 0; 2, -2, 0]; % avg position (reported)
n_tdoa = size(x_tdoa,2);

C_pos_full = .1*eye(2*n_tdoa); % position covar; all are IID
L = chol(C_pos_full,'lower');

% Generate a random set of TDOA positions
beta = x_tdoa + reshape(L*randn(2*n_tdoa,1),2,n_tdoa);

%% Generate Measurements
x_tgt = [6; 3];

zeta = tdoa.measurement(x_tdoa, x_tgt, n_tdoa);
zeta_unc = tdoa.measurement(beta, x_tgt, n_tdoa);

fprintf('Measurements from sensors 1-3 (w.r.t sensor 0):\n');
fprintf('Nominal Positions: %.2f m, %.2f m\n',zeta);
fprintf('Random Positions:  %.2f m, %.2f m\n',zeta_unc);

%% Plot Scenario
fig = figure;
scatter(x_tdoa(1,:),x_tdoa(2,:),'s','filled','DisplayName','Sensors (nominal positions)')
hold on;
scatter(beta(1,:), beta(2,:),'o','filled','DisplayName','Sensors (true positions)')
scatter(x_tgt(1), x_tgt(2),'^','filled','DisplayName','Target')
grid on;

% Draw the Isochrones
for idx=1:n_tdoa-1
    xy_iso = tdoa.drawIsochrone(x_tdoa(:,end),x_tdoa(:,idx),zeta_unc(idx),101,8);
    xy_iso_unc = tdoa.drawIsochrone(beta(:,end),beta(:,idx),zeta_unc(idx),101,8);

    hdl=plot(xy_iso(1,:),xy_iso(2,:),'-','DisplayName','Isochrone (nominal positions)');
    set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1);
    hdl_unc=plot(xy_iso_unc(1,:),xy_iso_unc(2,:),'--','DisplayName','Isochrone (true positions)');

    if idx>1
        utils.excludeFromLegend([hdl,hdl_unc]);
    end
end

legend('Location','SouthEast');
xlim([-1 8]);
ylim([-3 4]);

utils.setPlotStyle(gca,'widescreen');
