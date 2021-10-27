% Define Sensor Positions
baseline = [10e3,30e3];

nSensors = 3;
x_sensor1 = [-baseline(1), 0, baseline(1);0 0 0]/2;
x_sensor2 = [-baseline(2), 0, baseline(2);0 0 0]/2;

% Define Sensor Performance
timingError = [1e-7 1e-6];
% Ctoa = timingError^2*ones(nSensors,1);
Ctdoa1 = timingError(1)^2 * (1 + eye(nSensors-1));
Ctdoa2 = timingError(2)^2 * (1 + eye(nSensors-1));

% Define source positions
r_min = 0;
r_max = 1000e3;
M = 1001;
yvec = linspace(r_min,r_max,M);
x_source = [zeros(size(yvec));yvec];

% Compute CRLB
warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
crlb1 = tdoa.computeCRLB(x_sensor1,x_source,Ctdoa1); % Ndim x Ndim x M^2
crlb2 = tdoa.computeCRLB(x_sensor2,x_source,Ctdoa2); % Ndim x Ndim x M^2
cep1 = utils.computeCEP50(crlb1);
cep2 = utils.computeCEP50(crlb2);
warning('on','MATLAB:nearlySingularMatrix'); % Reactivate the singular matrix warning

figure;
loglog(yvec/1e3,cep1/1e3,yvec/1e3,cep2/1e3);
legend('10 km Baseline, 100 ns Timing Error','30 km Baseline, 1 \mus Timing Error');
xlabel('Distance [km]');ylabel('CEP_{50} [km]');
grid on;
