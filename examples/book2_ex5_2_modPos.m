% Executes a modified version of Example 5.2 from Practical Geolocation 
% for Electronic Warfare with MATLAB, used in Video 5.2.
%
% Nicholas O'Donoughue
% 8 February 2022

%% Set up sensors
alt1 = 1e3;
x_tdoa = [-15e3, -5e3, 5e3, 15e3; 0, 0, 0, 0; alt1, alt1, alt1, alt1];
[~,n_tdoa] = size(x_tdoa);

% Define target position and initial estimate
tgt_alt = 100; % known target altitude
x_tgt = [-10e3; 40e3; tgt_alt];
x_init = [0;10e3;alt1];

% Sensor Accuracy
time_err = 1e-7;
Croa = (utils.constants.c*time_err)^2 * eye(n_tdoa);
L = chol(Croa,'lower');

% Measurement and Noise
z = tdoa.measurement(x_tdoa, x_tgt, []);
noise = L*randn(n_tdoa,1); % generate sensor-level noise
noise_z = utils.resampleNoise(noise, []); % generate measurement-level noise
zeta = z + noise_z;

%% Solve for target position
[x_gd, x_gd_full] = tdoa.gdSoln(x_tdoa, zeta, Croa, x_init);

[a, ~] = utils.constraints.fixedAlt(tgt_alt, 'flat');
[x_gd_alt, x_gd_alt_full] = tdoa.gdSolnFixed(x_tdoa, zeta, Croa, x_init, a);

fprintf('Unconstrained Solution: %.2f km E, %.2f km N, %.2f km U\n',...
        x_gd(1)/1e3, x_gd(2)/1e3, x_gd(3)/1e3);

fprintf('Constrained Solution:   %.2f km E, %.2f km N, %.2f km U\n',...
        x_gd_alt(1)/1e3, x_gd_alt(2)/1e3, x_gd_alt(3)/1e3);

%% Plot the scenario
fig1 = figure;
stem3(x_tdoa(1,:),x_tdoa(2,:),x_tdoa(3,:),'k+','DisplayName','Sensors')
hold on;
stem3(x_tgt(1), x_tgt(2), x_tgt(3),'^','DisplayName','Target');
set(gca,'ColorOrderIndex',1);
grid on;

% Draw the Isochrones at alt=0
for tx_id = 1:n_tdoa-1
    iso = tdoa.drawIsochrone(x_tdoa(1:2,end),x_tdoa(1:2,tx_id),zeta(tx_id),101,40e3);
    hdl=plot(iso(1,:),iso(2,:),'--k','DisplayName','Isochrones');
    if tx_id > 1
        utils.excludeFromLegend(hdl);
    end
end

% Draw the GD
hdl=plot3(x_gd_full(1,:), x_gd_full(2,:), max(0,x_gd_full(3,:)),'-.');
utils.excludeFromLegend(hdl);
plot3(x_gd(1), x_gd(2), max(0,x_gd(3)),'-.s','Color',hdl.Color,'DisplayName','GD (unconstrained)');

hdl=plot3(x_gd_alt_full(1,:), x_gd_alt_full(2,:), max(0,x_gd_alt_full(3,:)),'-.');
utils.excludeFromLegend(hdl);
plot3(x_gd_alt(1), x_gd_alt(2), max(0,x_gd_alt(3)),'-.s','Color',hdl.Color,'DisplayName','GD (constrained');

xlim([-20e3,20e3]);
ylim([0e3,50e3]);
% zlim([0,12e3]);
legend();

view(-45,10);
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');


%% Retry with better elevation support
alt2=.5*alt1;
alt3=0*alt1;
x_tdoa2 = [-15e3, -5e3, 5e3, 15e3; 0, 5e3, 0, 5e3; alt2, alt1, alt3, alt2];

% Measurement and Noise
z2 = tdoa.measurement(x_tdoa2, x_tgt, []);
zeta2 = z2 + noise_z;

% Solve for target position
[x_gd2, x_gd2_full] = tdoa.gdSoln(x_tdoa2, zeta2, Croa, x_init);

% Plot old and new results
fig2 = figure;
stem3(x_tdoa(1,:),x_tdoa(2,:),x_tdoa(3,:),'+','DisplayName','Sensors')
hold on;
stem3(x_tdoa2(1,:),x_tdoa2(2,:),x_tdoa2(3,:),'+','DisplayName','Sensors (alt. config)')
stem3(x_tgt(1), x_tgt(2), x_tgt(3),'^','DisplayName','Target');
set(gca,'ColorOrderIndex',1);
grid on;

% Draw the Isochrones at alt=0
for tx_id = 1:n_tdoa-1
    iso = tdoa.drawIsochrone(x_tdoa2(1:2,end),x_tdoa2(1:2,tx_id),zeta2(tx_id),101,40e3);
    hdl=plot(iso(1,:),iso(2,:),'--k','DisplayName','Isochrones');
    if tx_id > 1
        utils.excludeFromLegend(hdl);
    end
end


% Draw the GD
hdl=plot3(x_gd_full(1,:), x_gd_full(2,:), max(0,x_gd_full(3,:)),'-.');
utils.excludeFromLegend(hdl);
plot3(x_gd(1), x_gd(2), max(0,x_gd(3)),'-.s','Color',hdl.Color,'DisplayName','GD (unconstrained)');

hdl=plot3(x_gd_alt_full(1,:), x_gd_alt_full(2,:), max(0,x_gd_alt_full(3,:)),'-.');
utils.excludeFromLegend(hdl);
plot3(x_gd_alt(1), x_gd_alt(2), max(0,x_gd_alt(3)),'-.s','Color',hdl.Color,'DisplayName','GD (constrained');

hdl=plot3(x_gd2_full(1,:), x_gd2_full(2,:), max(0,x_gd2_full(3,:)),'-.');
utils.excludeFromLegend(hdl);
plot3(x_gd2(1), x_gd2(2), max(0,x_gd2(3)),'-.s','Color',hdl.Color,'DisplayName','GD (alt. config)');

xlim([-20e3,20e3]);
ylim([0e3,50e3]);
% zlim([0,12e3]);
legend();

view(-45,10);
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');

fprintf('Unconstrained Solution (alt config):   %.2f km E, %.2f km N, %.2f km U\n',...
        x_gd2(1)/1e3, x_gd2(2)/1e3, x_gd2(3)/1e3);
%% Collect Figure Handles for Export
figs = [fig1, fig2];