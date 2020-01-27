% This script generates example data for homework Example 8.1 and
% Problems 8.5 and 8.6.
%
% INPUTS
%   none
%
% OUTPUTS
%   examples/ex8_1.mat      Sample data for Example 8.1
%   hw/problem8_5.mat       Sample data for Problem 8.5
%   hw/problem8_6.mat       Sample data for Problem 8.6
%
% Nicholas O'Donoughue
% 1 July 2019

% Initialize PTT Radio Parameters
Pt = 10;
Gt = 0;
Lt = 3;
B = 50e3;
ht = 6;
f0 = 750e6;

% Initialize receiver parameters
Gr = 0;
Lr = 2;
NF = 4;
Br = 1e6;
hr = 1e3;
N = 25;  % Array elements
M = 100; % Time samples
d_lam = .5;
v = array.make_steering_vector(d_lam,N);

% Initialize source positions
D = 3; % Sources
rng = [45 50 55]*1e3;
theta = [-30, 20, 24];
psi = theta*pi/180;

% Compute Received Power at each element
Lprop = prop.pathLoss(rng,f0,ht,hr,false,[]);
Pr = (10*log10(Pt) + Gt - Lt) - Lprop + Gr - Lr;
N_db = 10*log10(utils.constants.kT*Br)+NF;
SNR = Pr - N_db;
SNR_lin = 10.^(SNR/10); % D x 1
Pr_lin = 10.^(Pr/10);
N_lin = 10.^(N_db/10);

% Set up received data vector
V = v(psi);  % N x D
y = V*(sqrt(Pr_lin(:)/2).*(randn(D,M)+1i*(randn(D,M))));

% Add noise
n = sqrt(N_lin/2)*(randn(N,M)+1i*randn(N,M));
x = (y+n); % N x M

% Compute Array Factor
u_scan = linspace(-.8,.8,1001);
psi_scan = asin(u_scan);
th_scan = psi_scan*180/pi;
AF = array.compute_array_factor(v,x(:,1),psi_scan);
figure;
plot(th_scan,abs(AF));

save examples/ex8_1.mat x D N M d_lam;
% % Compute AoA Estimates
% [P,psi_vec] = array.beamscan(x,v,pi/2,1001);
% P_mvdr = array.beamscan_mvdr(x,v,pi/2,1001);
% 
% figure;
% th_vec = 180*psi_vec/pi;
% plot(th_vec,10*log10(abs(P)));
% hold on;
% set(gca,'ColorOrderIndex',3);
% plot(th_vec,10*log10(abs(P_mvdr)));
% legend('Beamscan','MVDR');
% xlabel('$\theta$ [deg]');
% ylabel('P [dB]');

%% Problem 8.5

% Initialize PTT Radio Parameters
Pt = 10;
Gt = 0;
Lt = 3;
B = 50e3;
ht = 6;
f0 = 750e6;

% Initialize receiver parameters
Gr = 0;
Lr = 2;
NF = 6;
Br = 1e6;
hr = 1e3;
N = 25;  % Array elements
M = 100; % Time samples
d_lam = .5;
v = array.make_steering_vector(d_lam,N);

% Initialize source positions
D = 3; % Sources
rng = [85 95 100]*1e3;
theta = [-10, 15, -5];
psi = theta*pi/180;

% Compute Received Power at each element
Lprop = prop.pathLoss(rng,f0,ht,hr,false,[]);
Pr = (10*log10(Pt) + Gt - Lt) - Lprop + Gr - Lr;
N_db = 10*log10(utils.constants.kT*Br)+NF;
SNR = Pr - N_db;
SNR_lin = 10.^(SNR/10); % D x 1
Pr_lin = 10.^(Pr/10);
N_lin = 10.^(N_db/10);

% Set up received data vector
V = v(psi);  % N x D
y = V*(sqrt(Pr_lin(:)/2).*(randn(D,M)+1i*(randn(D,M))));

% Add noise
n = sqrt(N_lin/2)*(randn(N,M)+1i*randn(N,M));
x = (y+n); % N x M

% Compute Array Factor
u_scan = linspace(-.8,.8,1001);
psi_scan = asin(u_scan);
th_scan = psi_scan*180/pi;
AF = array.compute_array_factor(v,x(:,1),psi_scan);
figure;
plot(th_scan,abs(AF));

save hw/problem8_5.mat x D N M d_lam;
% % Compute AoA Estimates
% [P,psi_vec] = array.beamscan(x,v,pi/2,1001);
% P_mvdr = array.beamscan_mvdr(x,v,pi/2,1001);
% 
% figure;
% th_vec = 180*psi_vec/pi;
% plot(th_vec,10*log10(abs(P)));
% hold on;
% set(gca,'ColorOrderIndex',3);
% plot(th_vec,10*log10(abs(P_mvdr)));
% legend('Beamscan','MVDR');
% xlabel('$\theta$ [deg]');
% ylabel('P [dB]');

%% Problem 8.6

% Initialize PTT Radio Parameters
Pt = 10;
Gt = 0;
Lt = 3;
B = 50e3;
ht = 6;
f0 = 750e6;

% Initialize receiver parameters
Gr = 0;
Lr = 5;
NF = 6;
Br = 1e6;
hr = 1e3;
N = 25;  % Array elements
M = 100; % Time samples
d_lam = .5;
v = array.make_steering_vector(d_lam,N);

% Initialize source positions
D = 3; % Sources
rng = [250 260 300]*1e3;
theta = [-10, 15, -5];
psi = theta*pi/180;

% Compute Received Power at each element
Lprop = prop.pathLoss(rng,f0,ht,hr,false,[]);
Pr = (10*log10(Pt) + Gt - Lt) - Lprop + Gr - Lr;
N_db = 10*log10(utils.constants.kT*Br)+NF;
SNR = Pr - N_db;
SNR_lin = 10.^(SNR/10); % D x 1
Pr_lin = 10.^(Pr/10);
N_lin = 10.^(N_db/10);

% Set up received data vector
V = v(psi);  % N x D
y = V*(sqrt(Pr_lin(:)/2).*(randn(D,M)+1i*(randn(D,M))));

% Add noise
n = sqrt(N_lin/2)*(randn(N,M)+1i*randn(N,M));
x = (y+n); % N x M

% Compute Array Factor
u_scan = linspace(-.8,.8,1001);
psi_scan = asin(u_scan);
th_scan = psi_scan*180/pi;
AF = array.compute_array_factor(v,x(:,1),psi_scan);
figure;
plot(th_scan,abs(AF));

save hw/problem8_6.mat x D N M d_lam;
% % Compute AoA Estimates
% [P,psi_vec] = array.beamscan(x,v,pi/2,1001);
% P_mvdr = array.beamscan_mvdr(x,v,pi/2,1001);
% 
% figure;
% th_vec = 180*psi_vec/pi;
% plot(th_vec,10*log10(abs(P)));
% hold on;
% set(gca,'ColorOrderIndex',3);
% plot(th_vec,10*log10(abs(P_mvdr)));
% legend('Beamscan','MVDR');
% xlabel('$\theta$ [deg]');
% ylabel('P [dB]');
