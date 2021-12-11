function fig=ex3_1()
% fig=ex3_1()
%
% Executes Example 3.1.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle to generated graphic
%
% Nicholas O'Donoughue
% 1 July 2019

% Check for existence of Statistics & Machine Learning Toolbox
use_stat_toolbox = license('test','Statistics_Toolbox');
   % If TRUE, then built-in functions will be used.
   % If FALSE, then custom-builts replacements in the utils namespace will
   % be used.
   
%  At each range, compute the received
%  power level.  Conduct a monte carlo trial
%  at each power level to compute PD and compare
%  to the desired threshold.

ht = 100;
hr = 2;
Rvec=100e3:1e3:1000e3;
Rvec_coarse = 100e3:50e3:1000e3;

f0 = 100e6;

% Compute Losses and Fresnel Zone
Lprop = prop.pathLoss(Rvec,f0,ht,hr,false);
Lprop_coarse = prop.pathLoss(Rvec_coarse,f0,ht,hr,false);

% Noise Power
B = 200e3; % channel bandwidth [Hz]
NF = 5; % noise figure [dB]
N0 = utils.constants.boltzmann*utils.constants.T0*10^(NF/10);
N = 10*log10(N0*B);

% Signal Power
EIRP = 47; % dBW
Gr = 0; % Receive antenna gain
Lr = 0;

% Received Power and SNR
Pr = EIRP-Lprop+Gr-Lr;
Pr_coarse = EIRP-Lprop_coarse+Gr-Lr;

% Desired Detection Performance
M = 10; % Number of samples
Pfa = 1e-6;
probD_desired=.5;

% Find expected max det range, for plotting
R_soln = detector.squareLawMaxRange(Pfa,probD_desired,M,f0,ht,hr,EIRP+Gr-Lr-N,false,[]);

% Monte Carlo Trial
nMC = 1e4;

% Convert noise and signal power to linear units
Nlin = 10.^(N/10);
Pr_lin = 10.^(Pr_coarse/10);
xi = Pr-N;
xi_lin = 10.^(xi/10);

% Theoretical Result
if use_stat_toolbox
    eta = chi2inv(1-Pfa,2*M);
    probD_theo = 1-ncx2cdf(eta,2*M,2*M.*xi_lin);
else
    eta = utils.chi2inv(1-Pfa,2*M);
    probD_theo = 1-utils.ncx2cdf(eta,2*M,2*M.*xi_lin);
end

% Generate noise and signal vectors
n = sqrt(Nlin/2) * (randn(M,nMC) + 1i*randn(M,nMC)); % Complex Gaussian with variance of N
start_phase = 2*pi*rand(1,nMC);
pulse = exp(1i*2*pi*(0:M-1)/10); % Unit Power
s = pulse(:).*exp(1i*start_phase);

% Compute Sufficient Statistic
probD = zeros(size(Pr_coarse));
for rngidx = 1:numel(Pr_coarse)
    % Scale signal power
    thisS = sqrt(Pr_lin(rngidx))*s;

    % Run Energy Detector
    detResult = detector.squareLaw(n+thisS,Nlin/2,Pfa);
    
    % Count detections
    probD(rngidx) = sum(detResult)/nMC;
end

fig=figure;
colors = get(0,'DefaultAxesColorOrder');
plot(Rvec/1e3,probD_theo,'Color',colors(2,:),'DisplayName','Theoretical');
hold on;
plot(Rvec_coarse/1e3,probD,'^','Color',colors(2,:),'DisplayName','Monte Carlo');
plot(R_soln/1e3*[1 1],[0 1],':','Color',colors(2,:),'DisplayName','Max Range');

plot(Rvec_coarse/1e3,probD_desired*ones(size(Rvec_coarse)),'k--','DisplayName','P_D=0.5');

xlabel('Range [km]');
ylabel('$P_D$');
%set(gca,'xscale','log');
ylim([0,1]);
legend('Location','SouthWest');

utils.setPlotStyle(gca,{'widescreen','tight'});
