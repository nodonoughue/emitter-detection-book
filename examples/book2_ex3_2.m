function figs = book2_ex3_2()
% figs=book2_ex3_2()
%
% Executes Example 3.2 from Practical Geolocation for Electronic Warfare
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
baseline = 10e3;
thSensors = linspace(0,2*pi,n_sensors) +pi/2; % add an extra sample, will be ignored
x_sensor = [0, cos(thSensors(1:end-1)); 0, sin(thSensors(1:end-1))] * baseline;
x_sensor = cat(1, x_sensor, 10e3*ones(1,n_sensors));

%% Define Covariance Matrix
time_err = 100e-9; % 100 ns
cov_full = time_err^2 * eye(n_sensors);

%% Plot Geometry
fig1=figure;
plot(x_sensor(1,:), x_sensor(2,:), 'o', 'DisplayName','Sensors');
legend('Location','NorthWest');

%% Sensor Pairs
% Cell array of the different configurations, to be iterated over
ref_set = {1, 'full', [1 3 3;2 4 2]};

%% Compute CRLB
figs = [fig1, zeros(size(ref_set))];

% Define search grid (targets up to 200 km away)
xx_vec = linspace(-100e3, 100e3, 101);
[xx,yy] = meshgrid(xx_vec);
alt=5e3;
x_source = [xx(:), yy(:), alt*ones(numel(xx),1)]';
%x_source = [xx(:), yy(:)]';

levels = [.01,1,5,10, 25, 50, 100, 200];

for idx_set = 1:numel(ref_set)
    this_ref = ref_set{idx_set};

    this_crlb = tdoa.computeCRLB(x_sensor, x_source, cov_full, this_ref);
        % N x N x 3
    
    this_cep = sqrt(this_crlb(1,1,:) + this_crlb(2,2,:) + this_crlb(3,3,:));
    %this_cep = utils.computeCEP50(this_crlb);
    
    % Plot this result
    figs(idx_set+1) = figure;
    [C,h] = contourf(xx_vec,xx_vec,reshape(this_cep/1e3, size(xx)),levels);
    clabel(C,h,'Color','white');
    set(gca,'ydir','normal');
    set(gca,'ColorScale','log');
    colormap(utils.viridis);
    caxis([levels(1),levels(end)]);
    colorbar
    hold on;
    plot(x_sensor(1,:), x_sensor(2,:), 'wo', 'DisplayName','Sensors');
end

%% Repeat with higher error for sensor one
cov_full(1,1) = 10*cov_full(1,1);
ref_set = {1,'full'};

figs2 = zeros(size(ref_set));

for idx_set = 1:numel(ref_set)
    this_ref = ref_set{idx_set};

    this_crlb = tdoa.computeCRLB(x_sensor, x_source, cov_full, this_ref);
        % N x N x 3
    
%    this_cep = sqrt(this_crlb(1,1,:) + this_crlb(2,2,:) + this_crlb(3,3,:));
    this_cep = utils.computeCEP50(this_crlb);
    
    % Plot this result
    figs2(idx_set) = figure;
    [C,h] = contourf(xx_vec,xx_vec,reshape(this_cep/1e3, size(xx)),levels);
    clabel(C,h,'Color','white');
    set(gca,'ydir','normal');
    set(gca,'ColorScale','log');
    colormap(utils.viridis);
    colorbar
    caxis([levels(1),levels(end)]);
    hold on;
    plot(x_sensor(1,:), x_sensor(2,:), 'wo', 'DisplayName','Sensors');
end