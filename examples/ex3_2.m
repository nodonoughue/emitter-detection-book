function fig=ex3_2()
% fig=ex3_2()
%
% Executes Example 3.2.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle to generated graphic
%
% Nicholas O'Donoughue
% 1 July 2019

%  At each range, compute the received
%  power level.  Conduct a monte carlo trial
%  at each power level to compute PD and compare
%  to the desired threshold.

% Check for existence of Statistics & Machine Learning Toolbox
use_stat_toolbox = license('test','Statistics_Toolbox');
   % If TRUE, then built-in functions will be used.
   % If FALSE, then custom-builts replacements in the utils namespace will
   % be used.
   
% Tranmsit Parameters
f0 = 35e9;
Pt = 10*log10(.035); % 35 mW / -14.56 dBW
Gt = 34; % dBi
Gr = 10;
Lt=0;
Lr = 2;
SLL = 25; % Sidelobe level (below Gt)
BW = 40e6;
NF = 4;
ht = 500;
hr = 500;

% Compute Noise Level and MDS
N0 = NF + 10*log10(utils.constants.boltzmann * utils.constants.T0);
N = N0 + 10*log10(BW);

% Compute Path Loss
Rvec = 100:25:100e3;
Rvec_coarse = [100:200:900,1e3:1e3:10e3,20e3:10e3:100e3];

alt = 1000;

Lprop = prop.pathLoss(Rvec,f0,alt,alt,true,[]);
Lprop_coarse = prop.pathLoss(Rvec_coarse,f0,alt,alt,true,[]);

% Compute Received Power
Pr = Pt + Gt - Lt+Gr-Lr-Lprop;
Pr_sl = Pr - SLL;
Pr_coarse = Pt + Gt - Lt+Gr-Lr-Lprop_coarse;
Pr_sl_coarse = Pr_coarse - SLL;

% Desired Detection Performance
M = 10; % Number of samples
Pfa = 1e-6;
probD_desired=.5;
SNRmin = detector.squareLawMinSNR(Pfa,probD_desired,M);
Pr_min = N + SNRmin;

% Find expected max det range, for plotting
R_soln = detector.squareLawMaxRange(Pfa,probD_desired,M,f0,ht,hr,Pt+Gt-Lt+Gr-Lr-N,true,[]);
R_sl_soln = detector.squareLawMaxRange(Pfa,probD_desired,M,f0,ht,hr,Pt+Gt-Lt+Gr-Lr-N-SLL,true,[]);

% Monte Carlo Trial
nMC = 1e4;


% Convert noise and signal power to linear units
Nlin = 10.^(N/10);
Pr_lin = 10.^(Pr_coarse/10);
Pr_sl_lin = 10.^(Pr_sl_coarse/10);

xi = Pr-N;
xi_lin = 10.^(xi/10);
xi_sl = Pr_sl-N;
xi_sl_lin = 10.^(xi_sl/10);

% Theoretical Result
if use_stat_toolbox
    eta = chi2inv(1-Pfa,2*M);
    probD_theo = 1-ncx2cdf(eta,2*M,2*M.*xi_lin);
    probD_sl_theo = 1-ncx2cdf(eta,2*M,2*M.*xi_sl_lin);
else
    eta = utils.chi2inv(1-Pfa,2*M);
    probD_theo = 1-utils.ncx2cdf(eta,2*M,2*M.*xi_lin);
    probD_sl_theo = 1-utils.ncx2cdf(eta,2*M,2*M.*xi_sl_lin);
end

% Generate noise and signal vectors
n = sqrt(Nlin/2) * (randn(M,nMC) + 1i*randn(M,nMC)); % Complex Gaussian with variance of N
start_phase = 2*pi*rand(1,nMC);
pulse = exp(1i*2*pi*(0:M-1)/10); % Unit Power
s = pulse(:).*exp(1i*start_phase);

% Compute Sufficient Statistic
probD = zeros(size(Pr_coarse));
probD_sl = zeros(size(Pr_coarse));
for rngidx = 1:numel(Pr_coarse)
    % Scale signal power
    thisS = sqrt(Pr_lin(rngidx))*s;
    thisS_sl = sqrt(Pr_sl_lin(rngidx))*s;
    
    % Run Energy Detector
    detResult = detector.squareLaw(n+thisS,Nlin/2,Pfa);
    detResult_sl = detector.squareLaw(n+thisS_sl,Nlin/2,Pfa);
    
    % Count detections
    probD(rngidx) = sum(detResult)/nMC;
    probD_sl(rngidx) = sum(detResult_sl)/nMC;
end


fig=figure;
colors = get(0,'DefaultAxesColorOrder');
plot(Rvec/1e3,probD_theo,'Color',colors(2,:),'DisplayName','Theoretical (ML)');
hold on;
plot(Rvec_coarse/1e3,probD,'^','Color',colors(2,:),'DisplayName','Monte Carlo (ML)');
plot(R_soln/1e3*[1 1],[0 1],':','Color',colors(2,:),'DisplayName','Max Range (ML)');

plot(Rvec/1e3,probD_sl_theo,'Color',colors(3,:),'DisplayName','Theoretical (SL)');
plot(Rvec_coarse/1e3,probD_sl,'^','Color',colors(3,:),'DisplayName','Monte Carlo (SL)');
plot(R_sl_soln/1e3*[1 1],[0 1],':','Color',colors(3,:),'DisplayName','Max Range (SL)');

plot(Rvec_coarse/1e3,probD_desired*ones(size(Rvec_coarse)),'k--','DisplayName','P_D=0.5');

xlabel('Range [km]');
ylabel('$P_D$');
set(gca,'xscale','log');
ylim([.01,1]);
legend('Location','South');

utils.setPlotStyle(gca,{'widescreen','tight'});
