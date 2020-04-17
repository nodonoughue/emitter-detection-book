function fig=ex10_2()
% fig=ex10_2()
%
% Executes Example 10.2 and generates one figure
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
x_sensor = 10*[-1 0;1 0]';

% Define measurement accuracy
sigma_psi = 2.5*pi/180;
C_psi = sigma_psi^2*eye(size(x_sensor,2)); % N x N identity matrix

% Find maximum cross-range position at 100 km downrange
x_source = [-100:100;100*ones(1,201)];
CRLB = triang.computeCRLB(x_sensor*1e3,x_source*1e3,C_psi);
cep50 = utils.computeCEP50(CRLB);

good_points = cep50 <= 25e3;
max_cross_range = max(abs(x_source(1,good_points)))*1e3;

x_max = 100;
x_vec = -x_max:.25:x_max;
[X,Y] = ndgrid(x_vec);
x0 = [X(:) Y(:)]';



% Compute CRLB
CRLB = triang.computeCRLB(x_sensor*1e3,x0*1e3,C_psi);
cep50 = reshape(utils.computeCEP50(CRLB),size(X)); % m

% Blank out y=0
nan_mask = abs(Y)<1e-6;
cep50(nan_mask)=NaN;

% Plot
fig=figure;
plot(x_sensor(1,:),x_sensor(2,:),'ko','DisplayName','AOA Sensors');
hold on;grid off;
contourLevels = [1 5 10 25];
[cp,hiso2]=contour(X,Y,cep50/1e3,contourLevels,'k');
utils.excludeFromLegend(hiso2);
clabel(cp,hiso2);
xlabel('x [km]');
ylabel('y [km]');
%title('Triangulation CRLB RMSE [km]');
legend('Location','NorthEast');

ylim([-100 100]);
utils.setPlotStyle(gca,{'equal'});
