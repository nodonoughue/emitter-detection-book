function figs = book2_ex5_4()
% fig=book2_ex5_4()
%
% Executes Example 5.4 from Practical Geolocation for Electronic Warfare
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
ref_lla = [25, -15, 0];
x_aoa_enu = [0, 50e3, 0;
             0, 0,    50e3;
            10, 10,  10]; % ENU coordinates
[x,y,z] = utils.enu2ecef(x_aoa_enu(1,:), x_aoa_enu(2,:), x_aoa_enu(3,:),...
                         ref_lla(1), ref_lla(2), ref_lla(3), 'deg','m');
x_aoa = [x(:), y(:), z(:)]'; % ECEF

n_aoa = size(x_aoa,2);

sat_lla = [27, -13, 575e3];
[x, y, z] = utils.lla2ecef(sat_lla(1), sat_lla(2), sat_lla(3), 'deg','m');
x_tgt = [x, y, z]'; % ECEF
[e, n, u] = utils.lla2enu(sat_lla(1), sat_lla(2), sat_lla(3), ...
                          ref_lla(1), ref_lla(2), ref_lla(3), 'deg', 'm');
x_tgt_enu = [e, n, u]';

%% Constraints
alt_low = 500e3;
alt_high = 600e3;
b = utils.constraints.boundedAlt(alt_low, alt_high, 'ellipse');
a = utils.constraints.fixedAlt(sat_lla(end),'ellipse');

%% Measurements
% Error
err_aoa = 3*pi/180; % 1 deg
C_aoa = err_aoa^2*eye(2*n_aoa);
L = chol(C_aoa,'lower');

z = triang.measurement(x_aoa, x_tgt);
n = L*randn(2*n_aoa,1);
zeta = z+n;

%% Solvers
[x,y,z] = utils.lla2ecef(ref_lla(1), ref_lla(2), 500e3);
x_init = [x, y, z]';

[x_gd, x_gd_full] = triang.gdSoln(x_aoa, zeta, C_aoa, x_init);
[x_gd_bnd, x_gd_bnd_full] = triang.gdSolnBounded(x_aoa, zeta, C_aoa, x_init, b);

x_ctr = x_init;
x_offset = [300e3, 300e3, 300e3];
x_step = [4e3, 4e3, 4e3];

[x_ml, ~, ~] = triang.mlSoln(x_aoa, zeta, C_aoa, x_ctr, x_offset, x_step);
[x_ml_bnd, ~, ~] = triang.mlSolnConstrained(x_aoa, zeta, C_aoa, x_ctr, x_offset, x_step, [], b);
[x_ml_fix, ~, ~] = triang.mlSolnConstrained(x_aoa, zeta, C_aoa, x_ctr, x_offset, x_step, a, [], 3e3);

%% Convert to LLA and Plot
[lat,lon,alt] = utils.ecef2lla(x_gd(1), x_gd(2), x_gd(3), 'm');
fprintf('Unconstrained Solution:  %.2f deg N, %.2f deg W, %.2f km\n', ...
        lat, abs(lon), alt/1e3);
fprintf('   Error: %.2f km\n', norm(x_gd-x_tgt)/1e3);
[lat,lon,alt] = utils.ecef2lla(x_gd_bnd(1), x_gd_bnd(2), x_gd_bnd(3), 'm');
fprintf('Constrained Solution:    %.2f deg N, %.2f deg W, %.2f km\n', ...
        lat, abs(lon), alt/1e3);
fprintf('  Error:  %.2f km\n', norm(x_gd_bnd-x_tgt)/1e3);

[lat,lon,alt] = utils.ecef2lla(x_ml(1), x_ml(2), x_ml(3), 'm');
fprintf('Unconstrained ML Solution:  %.2f deg N, %.2f deg W, %.2f km\n', ...
        lat, abs(lon), alt/1e3);
fprintf('   Error: %.2f km\n', norm(x_ml-x_tgt)/1e3);

[lat,lon,alt] = utils.ecef2lla(x_ml_bnd(1), x_ml_bnd(2), x_ml_bnd(3), 'm');
fprintf('Constrained ML Solution:  %.2f deg N, %.2f deg W, %.2f km\n', ...
        lat, abs(lon), alt/1e3);
fprintf('   Error: %.2f km\n', norm(x_ml_bnd-x_tgt)/1e3);

[lat,lon,alt] = utils.ecef2lla(x_ml_fix(1), x_ml_fix(2), x_ml_fix(3), 'm');
fprintf('Fixed Alt ML Solution:  %.2f deg N, %.2f deg W, %.2f km\n', ...
        lat, abs(lon), alt/1e3);
fprintf('   Error: %.2f km\n', norm(x_ml_fix-x_tgt)/1e3);

%% Plot in ENU
[e,n,u] = utils.ecef2enu(x_gd(1), x_gd(2), x_gd(3), ...
                         ref_lla(1), ref_lla(2), ref_lla(3), 'deg', 'm');
x_gd_enu = [e(:), n(:), u(:)]';

[e,n,u] = utils.ecef2enu(x_gd_full(1,:), x_gd_full(2,:), x_gd_full(3,:), ...
                         ref_lla(1), ref_lla(2), ref_lla(3), 'deg', 'm');
x_gd_full_enu = [e(:), n(:), u(:)]';

[e,n,u] = utils.ecef2enu(x_gd_bnd(1), x_gd_bnd(2), x_gd_bnd(3), ...
                         ref_lla(1), ref_lla(2), ref_lla(3), 'deg', 'm');
x_gd_bnd_enu = [e(:), n(:), u(:)]';

[e,n,u] = utils.ecef2enu(x_gd_bnd_full(1,:), x_gd_bnd_full(2,:), x_gd_bnd_full(3,:), ...
                         ref_lla(1), ref_lla(2), ref_lla(3), 'deg', 'm');
x_gd_bnd_full_enu = [e(:), n(:), u(:)]';

[e,n,u] = utils.ecef2enu(x_ml(1), x_ml(2), x_ml(3), ...
                         ref_lla(1), ref_lla(2), ref_lla(3), 'deg', 'm');
x_ml_enu = [e, n, u]';

[e,n,u] = utils.ecef2enu(x_ml_bnd(1), x_ml_bnd(2), x_ml_bnd(3), ...
                         ref_lla(1), ref_lla(2), ref_lla(3), 'deg', 'm');
x_ml_bnd_enu = [e, n, u]';

fig1=figure;
stem3(x_aoa_enu(1,:),x_aoa_enu(2,:),x_aoa_enu(3,:),'o','filled','DisplayName','Sensors')
hold on;
stem3(x_tgt_enu(1), x_tgt_enu(2), x_tgt_enu(3),'^','DisplayName','Target');
set(gca,'ColorOrderIndex',1);
grid on;

% Draw the GD
hdl=plot3(x_gd_full_enu(1,:), x_gd_full_enu(2,:), max(0,x_gd_full_enu(3,:)),'-.');
utils.excludeFromLegend(hdl);
plot3(x_gd_enu(1), x_gd_enu(2), max(0,x_gd_enu(3)),'-.o','Color',hdl.Color,'DisplayName','GD (unconstrained)');

hdl=plot3(x_gd_bnd_full_enu(1,:), x_gd_bnd_full_enu(2,:), max(0,x_gd_bnd_full_enu(3,:)),'-.');
utils.excludeFromLegend(hdl);
plot3(x_gd_bnd_enu(1), x_gd_bnd_enu(2), max(0,x_gd_bnd_enu(3)),'-.s','Color',hdl.Color,'DisplayName','GD (constrained');

stem3(x_ml(1), x_ml(2), x_ml(3),'v','DisplayName','ML (unconstrained)');
stem3(x_ml_bnd(1), x_ml_bnd(2), x_ml_bnd(3), '+','DisplayName','ML (constrained)');
legend();

view(-45,10);
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');

%% Collect Figure Handles for Export
figs = [fig1];