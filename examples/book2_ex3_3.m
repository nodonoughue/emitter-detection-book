function figs = book2_ex3_3()
% figs=book2_ex3_3()
%
% Executes Example 3.3 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   figs         array of figure handles
%
% Nicholas O'Donoughue
% 31 May 2021

%% Set up sensor coordinates
n_sensors = 4;
alt=10e3;
x_sensor = [-17.5e3, -12.5e3, 12.5e3, 17.5e3;
            zeros(1,n_sensors)];
%            alt*ones(1,n_sensors)];
vel = 100; % m/s
th_1 = 70 * pi/180;
th_2 = 110 * pi/180;
v_sensor = [cos(th_1)*ones(1,2), cos(th_2)*ones(1,2);
            sin(th_1)*ones(1,2), sin(th_2)*ones(1,2)];
%            zeros(1,n_sensors)]*vel;
      
%% Define Covariance Matrix
lam = 3e8/1e9;
time_err = 100e-9; % 100 ns
freq_err = 100; % 10 Hz
rng_err = utils.constants.c*time_err;
rng_rate_err = lam*freq_err;
cov_tdoa = rng_err^2 * eye(n_sensors);
cov_fdoa = rng_rate_err^2 * eye(n_sensors);
cov_full = blkdiag(cov_tdoa, cov_fdoa);

%% Plot Geometry
fig1=figure;
plot(x_sensor(1,:), x_sensor(2,:), 'o', 'DisplayName','Sensors');
hold on;
for idx_sensor=1:n_sensors
    utils.drawArrow(x_sensor(1,idx_sensor) + [0 v_sensor(1,idx_sensor)]*1e3, ...
                    x_sensor(2,idx_sensor) + [0 v_sensor(2,idx_sensor)]*1e3);
end
ylim([-10e3 10e3])
legend('Location','NorthWest');
utils.setPlotStyle(gca,'equal');

%% Sensor Pairs
% Cell array of the different configurations, to be iterated over
ref_set = {[1 3;2 4], [1 2 3;2 3 4], 'full'};

%% Compute CRLB
figs = [fig1, zeros(size(ref_set))];

% Define search grid (targets up to 200 km away)
xx_vec = linspace(-50e3, 50e3, 301);
[xx,yy] = meshgrid(xx_vec);
x_source = [xx(:), yy(:)]';%, alt*ones(numel(xx),1)]';

levels = [.1,.5,1,5,10,25];

for idx_set = 1:numel(ref_set)
    this_ref = ref_set{idx_set};

    % Repeat x_sensors for both the TDOA and FDOA sensor sets, because all
    % sensors do both.  Repeat the ref index, because the sensor pairs for
    % TDOA and FDOA are the same.
    this_crlb = hybrid.computeCRLB([],x_sensor, x_sensor, v_sensor, ...
                                   x_source, cov_full, this_ref, this_ref);
        % N x N x 3
    
    this_rmse = sqrt(this_crlb(1,1,:) + this_crlb(2,2,:));% + this_crlb(3,3,:));

    % Plot this result
    figs(idx_set+1) = figure;
    [C,h] = contourf(xx_vec,xx_vec,reshape(this_rmse/1e3, size(xx)),levels);
    clabel(C,h,'Color','white');
    set(gca,'ydir','normal');
    colormap(utils.viridis);
    set(gca,'ColorScale','log');
    colorbar
    hold on;
    plot(x_sensor(1,:), x_sensor(2,:), 'wo', 'DisplayName','Sensors');
    for idx_sensor=1:n_sensors
        utils.drawArrow(x_sensor(1,idx_sensor) + [0 v_sensor(1,idx_sensor)]*1e1, ...
                        x_sensor(2,idx_sensor) + [0 v_sensor(2,idx_sensor)]*1e1);
    end
end


