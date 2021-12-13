% Draw Figures - Chapter 3
%
% This script generates all of the figures that appear in
% Chapter 3 of the textbook.
%
% Nicholas O'Donoughue
% 1 July 2019

% Clear Figures
close all;

% Set up directory and filename for figures
dirNm = fullfile(pwd,'figures');
if ~exist(dirNm,'dir')
    mkdir(dirNm);
end
prefix = fullfile(dirNm,'fig3_');

% Initialize Plot Preference
utils.initPlotSettings;

% Initializes colorSet - Mx3 RGB vector for successive plot lines
colors = get(0,'DefaultAxesColorOrder');

% Reset the random number generator, to ensure reproducability
rng('default') ; 

% Add path to folder for examples from the textbook
addpath('examples');

% Check for existence of Statistics & Machine Learning Toolbox
use_stat_toolbox = license('test','Statistics_Toolbox');
   % If TRUE, then built-in functions will be used.
   % If FALSE, then custom-builts replacements in the utils namespace will
   % be used.
   
%% Figure 1, Spectral Content

% Frequency of Signal
f0 = 1;
BW = .4;

% Amplitudes
N0 = 1;
A = 5;

% Generate Frequency Content
N = 201;
f = 2*f0*linspace(-1,1,N);
n = N0*ones(1,N);
s = max(0,A*(1-2*abs(abs(f)-f0)/BW));

% Plot
fig1=figure;
set(gca,'ColorOrder',colors(1:3:end,:));hold on;
plot(f,n,'DisplayName','Noise');
plot(f,s,'DisplayName','Signal');

xlabel('$f$');
ylabel('$P(f)$');
legend('Location','NorthEast');

utils.setPlotStyle(gca,{'box only','widescreen','tight'});
utils.exportPlot(fig1,[prefix '1']);


%% Figure 2 - Plot with Spectrum

BWf = 2*BW;
H = zeros(1,N) + 1.2*A*(abs(abs(f)-f0)<=(BW/2));
n_filtered = n.*(H>0);

% Plot
fig2=figure;
set(gca,'ColorOrder',colors(1:3:end,:));
plot(f,n,'DisplayName','Noise');hold on;
plot(f,s,'DisplayName','Signal');
plot(f,H,'--','DisplayName','Filter');

xlabel('$f$');
ylabel('$P(f)$');
legend('Location','NorthEast');
ylim([0 ceil(1.2*max(H(:)))]);

utils.setPlotStyle(gca,{'box only','tight'});

utils.exportPlot(fig2,[prefix '2']);

%% Figure 3, CW Detection PFA vs. Threshold
M = [1 10 100];
eta_db = -10:.1:30;
eta_lin = 10.^(eta_db/10);

% Generate ND grid for easy calculation
[MM,ETA] = ndgrid(M,eta_lin);

% Compute PFA for each eta,N pair
if use_stat_toolbox
    Pfa = chi2cdf(ETA,2*MM,'upper');
else
    Pfa = utils.chi2cdf(ETA,2*MM,'upper');
end


% Plot
fig3=figure;
semilogy(eta_db,Pfa);
legendStr = arrayfun(@(x) sprintf('M=%d',x),M,'UniformOutput',false);
legend(legendStr,'Location','SouthWest');
xlabel('$\eta$ [dB]');
ylabel('$P_{FA}$');
grid on;
ylim([1e-6 1]);xlim([eta_db(1),eta_db(end)]);

utils.setPlotStyle(gca,{'widescreen','tight'});

utils.exportPlot(fig3,[prefix '3']);

%% Figure 4, PD vs. SNR for CW Detection
PFA = 1e-6;
M = [1 10 100 1000];
xi_db = -20:.1:20;
xi_lin = 10.^(xi_db/10);

% Use NDGrid to generate N,xi pairs
[MM,XI] = ndgrid(M,xi_lin);

% Compute threshold
if use_stat_toolbox
    eta = chi2inv(1-PFA,2*MM);
else
    eta = utils.chi2inv(1-PFA,2*MM);
end

% Compute Probability of Detection
lambda = 2.*XI;
if use_stat_toolbox
    PD = 1 - ncx2cdf(eta,2*MM,MM.*lambda);
else
    PD = 1 - utils.ncx2cdf(eta,2*MM,MM.*lambda);
end

% Plot
fig4 = figure;
plot(xi_db,PD);
legendStr = arrayfun(@(x) sprintf('M=%d',x),M,'UniformOutput',false);
legend(legendStr,'Location','NorthWest');
xlabel('$\xi$ [dB]');
ylabel('$P_D$');
grid on;

utils.setPlotStyle(gca,{'widescreen','tight'});

utils.exportPlot(fig4,[prefix '4']);

%% Figure 7, Atmospheric Loss Table
R=1e3; % set ref distance to 1 km
f = 1e9:50e6:100e9;

% Reference Atmosphere
%  -- Sea Level, 10 kft, 20 kft, 30 kft, 40 kft
alt_kft = [0,10,20,30,40]; % kft
legendCell = arrayfun(@(x) sprintf('%d kft',x),alt_kft,'UniformOutput',false);
% T = [15, -4.8, -24.6, -44.4, -56.6];
% P = [101325, 69680, 46560, 30090,18750];
% g = [7.5,2.01,0.34,.05,.01];

atmLoss = zeros(numel(f),numel(alt_kft));
for ii=1:numel(alt_kft)
    atmStruct = atm.standardAtmosphere(alt_kft(ii)*unitsratio('m','ft')*1e3);
    atmLoss(:,ii) = atm.calcAtmLoss(f,R,0,0,atmStruct);
end

fig7=figure;
semilogy(f/1e9,atmLoss(:,:));
legend(legendCell,'Location','NorthWest');
xlabel('Frequency [GHz]');
ylabel('Specific Attenuation [dB/km]');
grid on;

utils.setPlotStyle(gca,{'widescreen','tight'});

utils.exportPlot(fig7,[prefix '7']);


%% Figures 8, FM Reception Power vs. Range
% Figure 8 : SNR vs. range

% Set up RF environment
ht = 100;
hr = 2;
% Rvec = [1e3:1e3:9e3,10e3:10e3:90e3,100e3:100e3:1000e3];
Rvec=10e3:10e3:1000e3;
f0 = 100e6;

% Compute Losses and Fresnel Zone
Lfspl = prop.freeSpacePathLoss(Rvec,f0,false);
Ltworay = prop.twoRayPathLoss(Rvec,f0,ht,hr,false);
Lprop = prop.pathLoss(Rvec,f0,ht,hr,false);

% Noise Power
B = 2e5; % channel bandwidth [Hz]
NF = 5; % noise figure [dB]
N0 = utils.constants.boltzmann*utils.constants.T0*10^(NF/10);
N = 10*log10(N0*B);

% Signal Power
EIRP = 47; % dBW
Gr = 0; % Receive antenna gain
Lr = 0;

% Received Power and SNR
Pr = EIRP-Lprop+Gr-Lr;
SNRmin = 3.65;
Pr_min = N+SNRmin;

SNR0 = EIRP+Gr-Lr-N;
Rsoln = detector.squareLawMaxRange(1e-6,.5,10,f0,ht,hr,SNR0,false,[]);

fig8=figure;
plot(Rvec/1e3,Pr,'DisplayName','P_R');
hold on;
plot(Rvec/1e3,Pr_min*ones(size(Rvec)),':','DisplayName','MDS');
legend('Location','NorthEast');
xlabel('Range [km]');
ylabel('Received Power [dBW]');
grid on

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig8,[prefix '8']);

%% Figure 9, Monte Carlo Results
fig9 = ex3_1;
utils.exportPlot(fig9,[prefix '9']);

%% Figure 10, Monte Carlo Results
fig10= ex3_2;
utils.exportPlot(fig10,[prefix '10']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;