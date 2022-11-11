function fig = ex10_3()
% fig = ex10_3()
%
% Executes Example 10.3 and generates one figure.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle to generated graphic
%
% Nicholas O'Donoughue
% 1 July 2019

% Define sensor positions
x_sensor = 10*[-1 0;1 0;0 1]';
    
x_max = 100;
x_vec = -x_max:.25:x_max;
[X,Y] = ndgrid(x_vec);
x0 = [X(:) Y(:)]';

% Define measurement accuracy
sigma_psi = 2.5*pi/180;
C_psi = sigma_psi^2*eye(size(x_sensor,2)); % N x N identity matrix

% Compute CRLB
CRLB = triang.computeCRLB(x_sensor*1e3,x0*1e3,C_psi);
cep50 = reshape(utils.computeCEP50(CRLB),size(X)); % m

good_point = cep50 <= 25e3;
rng_val = sqrt(sum(abs(x0).^2,1)); % km

max_range = min(rng_val(~good_point)); % km
fprintf('Max range that satisfies CEP < 25 km regardless of AOA is: %.2f\n',max_range);


% Plot
fig=figure;
plot(x_sensor(1,:),x_sensor(2,:),'ko','DisplayName','AOA Sensors');
hold on;grid off;
contourLevels = [.1 .5 1 5 10 25];
[cp,hiso2]=contour(X,Y,cep50/1e3,contourLevels,'k');
utils.excludeFromLegend(hiso2);
clabel(cp,hiso2);
xlabel('x [km]');
ylabel('y [km]');
% title('Triangulation CRLB RMSE [km]');
legend('Location','NorthEast');

utils.setPlotStyle(gca,{'equal'});
