% Draw Figures - Chapter 2
%
% This script generates all of the figures that appear in
% Chapter 2 of the textbook.
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
prefix = fullfile(dirNm,'fig2_');

% Initialize Plot Preference
utils.initPlotSettings;

% Initializes colorSet - Mx3 RGB vector for successive plot lines
colors = get(0,'DefaultAxesColorOrder');

% Reset the random number generator, to ensure reproducability
rng('default') ; 

% Add path for folder with examples from the textbook
addpath('examples');

%% Figure 1, Likelihood Function

% Number of x axis points
N = 512;

% Centroids
n1 = .3;
n2 = .45;
n3 = .6;
s2 = .02; % variance of .5

% Generate PDFs
x = linspace(0,1,N);

% Use an exponential (unscaled) and a random factor to make them look a
% little less clean
f1 = exp(-(x-n1).^2/(2*s2));% + .05*randn(1,N);
f2 = exp(-(x-n2).^2/(2*s2));% + .05*randn(1,N);
f3 = exp(-(x-n3).^2/(2*s2));% + .05*randn(1,N);

% Scale PDFs
f1 = f1/norm(f1);
f2 = f2/norm(f2);
f3 = f3/norm(f3);

% Plot
fig1=figure;hold on;
set(gca,'ColorOrder',.1*ones(1,3));
plot(x,f1,'DisplayName','\theta=\theta_1');
plot(x,f2,'DisplayName','\theta=\theta_2');
plot(x,f3,'DisplayName','\theta=\theta_3');

xlabel('$z$');
ylabel('$f_\theta(z)$');
legend('Location','NorthEast');

utils.setPlotStyle(gca,{'box only','notick','widescreen','tight'});
utils.exportPlot(fig1,[prefix '1']);


%% Figure 2 - Likelihood Ratio Test

% Set up axes
ell = 0:.01:10;

% Set variance to 2 under both hypotheses, mean to 4 and 6 under H0 and H1,
% respectively
s2 = 2;
m0 = 2;
m1 = 6;

% PDF of ell(z) under H0 and H1, respectively
f0 = (2*pi*s2)^(-.5)*exp(-.5*(ell-m0).^2/s2);
f1 = (2*pi*s2)^(-.5)*exp(-.5*(ell-m1).^2/s2);

% Detection threshold
eta = 4.5;
mask = ell>=eta; % Region above the threshold
fa = f0(mask);
misseddet = f1(~mask);

% Plot the likelihood functions
fig2=figure;hold on;
set(gca,'ColorOrder',.1*ones(1,3));
plot(ell,f0,'DisplayName','H_0');
plot(ell,f1,'DisplayName','H_1');

% Add the threshold
plot(eta*[1 1],[0 .4],'k-.','LineWidth',.5,'DisplayName','\eta','Color',colors(2,:));
%text(eta+.1,.38,'$\eta$');

% Add the PFA/PD regions
area(ell(~mask) ,misseddet,'DisplayName','Missed Detection','FaceColor',colors(3,:),'FaceAlpha',.6,'LineStyle','none');
area(ell(mask),fa,'DisplayName','False Alarm','FaceColor',colors(4,:),'FaceAlpha',.6,'LineStyle','none');

% Add text overlay
annotation(fig2,'doublearrow',[0.433 0.516],...
    [0.885 0.885]);
text(eta+.5,.35,'Reduce $P_{FA}$');
text(eta-2.1,.35,'Reduce $P_{MD}$');

% Axis Labels
ylabel('$f_\theta(\ell)$');
xlabel('$\ell(z)$');

% Turn on the legend
legend('Location','NorthEast');

% Clean up the plot
utils.setPlotStyle(gca,{'box only','notick','widescreen','tight'});

% Output the files
utils.exportPlot(fig2,[prefix '2']);

%% Figure 3, Receiver Operating Characteristics
fig3 = ex2_2;
utils.exportPlot(fig3,[prefix '3']);

%% Figure 4, PD S-Curves

% Set up PFA and SNR vectors
pfa = [1e-6,1e-4,1e-2];
xi = 0:.1:10;           % dB Scale
xi_lin = 10.^(xi/10);   % Convert SNR to linear
legendStr = arrayfun(@(x) sprintf('P_{FA} = 10^{%d}',log10(x)),pfa,'UniformOutput',false);

% Use ndgrid to simplify the code
[PFA,XI] = ndgrid(pfa,xi_lin);

% Compute the PD according to the simplified
% erf/erfinv equation
PD = .5*(1-erf(erfinv(1-2*PFA)-10.^(XI/10)/sqrt(2)));

% Plot the ROC curve
fig4=figure;
plot(xi,PD);

% Axes Labels
ylabel('$P_D$');xlabel('$\xi$ [dB]');

% Legend
legend(legendStr,'Location','NorthWest');

% Align the axes
utils.setPlotStyle(gca,{'widescreen','tight'});


% Output the files
utils.exportPlot(fig4,[prefix '4']);

%% Figure 5, CFAR

% Noise Level
s2 = 1;
s2high = 10;

% Generate Noise
N = 1024; % number of samples
N1 = fix(N/2);   % start of high sample
N2 = fix(3*N/4); % end of high sample

n = sqrt(s2/2)*(randn(N,1)+1i*randn(N,1));
n(N1:N2) = sqrt(s2high/2)*(randn(N2-N1+1,1)+1i*randn(N2-N1+1,1));

% Variance Estimate
M = 64; % window size - one sided
Mg = 4; % Number of guard cells - one sided
mask = [ones(M,1); zeros(2*Mg+1,1); ones(M,1)]; % Initial CA-CFAR mask
mask = mask/(2*M);
s2est = conv(abs(n).^2,mask,'same');

% Threshold Equation --- need to replace
% eta = s2est*10;
eta = sum(toeplitz(s2est(1:100),s2est),1)/10;
% Plot
fig5=figure;hold on;
set(gca,'ColorOrder',.2*ones(1,3));
plot(20*log10(abs(n)),'DisplayName','Noise Level');
plot(10*log10(abs(eta)),'DisplayName','Threshold');

xlabel('Time/Frequency/Range');
ylabel('Power [dB]');

% Legend
legend('Location','NorthWest');

% Align the axes
utils.setPlotStyle(gca,{'box only','notick','widescreen','tight'});

% Output the files
utils.exportPlot(fig5,[prefix '5']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;