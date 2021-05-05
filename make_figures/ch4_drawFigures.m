% Draw Figures - Chapter 4
%
% This script generates all of the figures that appear in
% Chapter 4 of the textbook.
%
% Nicholas O'Donoughue
% 1 July 2019

% Clear Figures
close all;

% Set up directly/filename for figures
dirNm = fullfile(pwd,'figures');
dirNmComponents = fullfile(pwd,'Graphics','Components');
if ~exist(dirNm,'dir')
    mkdir(dirNm);
end
if ~exist(dirNmComponents,'dir')
    mkdir(dirNmComponents);
end

prefix = fullfile(dirNm,'fig4_');

% Initialize Plot Preference
utils.initPlotSettings;

% Initializes colorSet - Mx3 RGB vector for successive plot lines
colors = get(0,'DefaultAxesColorOrder');

% Reset the random number generator, to ensure reproducability
rng('default') ; 

% Add the examples folder to the path, for figures from
% examples
addpath('examples');

% Check for existence of Statistics & Machine Learning Toolbox
use_stat_toolbox = license('test','Statistics_Toolbox');
   % If TRUE, then built-in functions will be used.
   % If FALSE, then custom-builts replacements in the utils namespace will
   % be used.
   
%% Figure 1a - Alternating Sine Waves

% Sine wave
N = 1024; % Sample points
y = exp(1i*(pi/2+2*pi*(0:N-1)/N));

% String together multiple periods
code = [0 1 1 0 1];
T = numel(code)*N;
Y = zeros(1,T);
for ii=1:numel(code)
    Y((1:N)+N*(ii-1)) = y*exp(1i*pi*code(ii));
end

fig1a=figure;
plot(0:T-1,real(Y),[0,T-1],[0,0],'k-','LineWidth',.5);
hold on;
for ii=1:numel(code)
    text(N/2 + N*(ii-1),1.5,sprintf('%d',code(ii)));
end
for ii=0:numel(code)
    plot((N*ii)*[1 1],[-1,2],'k:');
end

% Annotation
annotation(fig1a,'doublearrow',[0.364087301587302 0.517857142857143],...
    [0.72563139329806 0.72663139329806]);
text(2.35*N,1.25,'$T_{chip}$');

utils.setPlotStyle(gca,{'widescreen','clean','tight'});
utils.exportPlot(fig1a,[prefix '1a']);

%% Figure 1b, Bandwidth

N = 16; % Number of samples per cycle
nCode = 128; % Length of transmit code in bits
y = exp(1i*(pi/2+2*pi*(0:N-1)/N));

nFull = nCode*N;
lo = exp(1i*2*pi*(0:nFull-1)*4/N);

nMC = 100;
XX = zeros(size(lo));
for ii=1:nMC
    code = randi(2,nCode)-1;
    phi1 = exp(1i*2*pi*rand(1));
    Y = zeros(1,nFull);
    for jj=1:nCode
        Y((1:N)+N*(jj-1)) = y.*exp(1i*pi*code(jj)).*phi1;
    end

    % Mix Y with a carrier at 8x the frequency
    X = Y.*lo;

    % Take the fourier transform
    XX = XX+abs(fftshift(fft((real(X)))));
end
XX = XX./max(abs(XX));

fig1b=figure;
plot(linspace(-1,1,nFull),10*log10(abs(XX)));
hold on;
% Plot top and -3 dB lines
plot([-1 1],[0 0],'k:');
plot([-1 1],[-3 -3],'k:');
plot([0 0],[-20 0],'k-');
ylim([-20,3]);

% Create textbox
text(-.4,-1.5,'3 dB');
text(.18,-4.5,'$B_s=1/T_{\mathrm{chip}}$');
text(.3,2,'$f_0$');

% Create doublearrows
annotation(fig1b,'doublearrow',[0.310515873015873 0.311507936507936],...
    [0.893179894179894 0.761904761904762]);
annotation(fig1b,'doublearrow',[0.77579365079365 0.846230158730158],...
    [0.766195767195767 0.767195767195767]);
annotation(fig1b,'doublearrow',[0.516865079365079 0.812499999999999],...
    [0.91963492063492 0.918871252204585]);
annotation(fig1b,'line',[0.739087301587302 0.814484126984127],...
    [0.706231040564374 0.756613756613757],'LineWidth',0.5);
utils.setPlotStyle(gca,{'clean','widescreen','tight'});
utils.exportPlot(fig1b,[prefix '1b']);

%% Figure 2a and 2b: Chip Rate and Bandwidth

% Generate the digital signals
Rdata = 4; % bits/sec
Rchip = 16; % bits/sec
nCode = Rchip/Rdata;
nData = 4;
nFull = nCode*nData;

dataBits = randi(2,1,nData)-1;
codeBits = randi(2,1,nCode)-1;

dataBitsFull = reshape(repmat(dataBits,nCode,1),1,nFull);
codeBitsFull = repmat(codeBits,1,nData);
outBitsFull = xor(dataBitsFull,codeBitsFull);

% Generate the signals
osf = 16; % Samples per cycle
y = exp(1i*(pi/2+2*pi*(0:osf-1)/osf));

% Construct the code signals
yData = reshape(y(:).*exp(1i*pi*dataBitsFull),1,nFull*osf);
yCode = reshape(y(:).*exp(1i*pi*codeBitsFull),1,nFull*osf);
yDSSS = reshape(y(:).*exp(1i*pi*outBitsFull),1,nFull*osf);

fig2a=figure;
% Start with the Signals at the origin
plot(1:nFull*osf,real(yData)+6);
hold on;
set(gca,'ColorOrderIndex',3)
plot(1:nCode*osf,real(yCode(1:nCode*osf))+3);
plot(1:nFull*osf,real(yDSSS));
hold on;
% Add the code and vertical lines
for i=1:nFull
    text(osf*(i-1)+osf/2,1.5,sprintf('%d',outBitsFull(i)));
    h=plot(osf*(i-1)*[1 1],[-1 2],'k:','LineWidth',.5);
    utils.excludeFromLegend(h);
end

for i=1:nCode
    text(osf*(i-1)+osf/2,4.5,sprintf('%d',codeBits(i)));
    h=plot(osf*(i-1)*[1 1],[2 5],'k:','LineWidth',.5);
    utils.excludeFromLegend(h);
end

for i=1:nData
    text(osf*nCode*(i-1)+osf*nCode/2,7.5,sprintf('%d',dataBits(i)));
    h=plot(osf*nCode*(i-1)*[1 1],[2 8],'k:','LineWidth',.5);
    utils.excludeFromLegend(h);
end
legend('Data Signal','Spreading Code','Encoded Signal','Location','East');
utils.setPlotStyle(gca,{'clean','widescreen','tight'});
utils.exportPlot(fig2a,[prefix '2a']);

nMC=1000;

nCode1 = 16;
nCode2 = 64;
nFull = 16*max(nCode1,nCode2); % at least 16 samples per period
N1 = nFull/nCode1;
N2 = nFull/nCode2;

% Generate the chips
y1 = (exp(1i*(pi/2+2*pi*(0:N1-1)/N1)));
y2 = (exp(1i*(pi/2+2*pi*(0:N2-1)/N2)));

YY1 = zeros(1,nFull);
YY2 = zeros(1,nFull);
for jj=1:nMC
    % Generate new codes
    c1 = randi(2,1,nCode1)-1;
    c2 = randi(2,1,nCode2)-1;
    
    phi1 = exp(1i*2*pi*rand(1));
    phi2 = exp(1i*2*pi*rand(1));
    % Construct the full signals
    Y1 = zeros(1,nFull);
    for ii=1:nCode1
        Y1((1:N1)+N1*(ii-1)) = y1*exp(1i*pi*c1(ii))*phi1;
    end
    Y2 = zeros(1,nFull);
    for ii=1:nCode2
        Y2((1:N2)+N2*(ii-1)) = y2*exp(1i*pi*c2(ii))*phi2;
    end
    
    % Take the fourier transform
    YY1 = YY1+abs(fftshift(fft(Y1)));
    YY2 = YY2+abs(fftshift(fft(Y2)));    
end
YY1 = abs(YY1)./max(abs(YY1));
YY2 = abs(YY2)./max(abs(YY2));

% Shift to account for central frequency of chip waveform
YY1 = circshift(YY1,-nCode1);
YY2 = circshift(YY2,-nCode2);
% Plot
fig2b = figure;
plot(linspace(0,1,nFull),10*log10(YY1));
hold on;
set(gca,'ColorOrderIndex',4);
plot(linspace(0,1,nFull),10*log10(YY2));
ylim([-20,2]);xlim([.25 .75]);
hold on;
legend('Data Signal','Encoded Signal');
% Plot top and -3 dB lines
h=plot([0 1],[-3 -3],'k:');
utils.excludeFromLegend(h);
text(.3,-2.5,'-3 dB');
text(.55,-2.25,'$B_s = 1/T_{chip}$');
text(.55,-4.5,'$B_d = 1/T_{sym}$');

% Annotation
annotation(fig2b,'line',[0.581597222222222 0.598090277777778],...
    [0.79220987654321 0.820987654320988]);
annotation(fig2b,'line',[0.605034722222222 0.529513888888889],...
    [0.733567901234568 0.78395061728395]);

utils.setPlotStyle(gca,{'clean','widescreen','tight'});
utils.exportPlot(fig2b,[prefix '2b']);
    
%% Figure 3, Spreading of SNR

snr_ref_db = 0;
bw_ref = 1e6;

bw = 10.^(3:.1:9);

spreadingLoss = 10*log10(bw./bw_ref);

snr_o_db = snr_ref_db - spreadingLoss;

fig3=figure;
semilogx(bw,snr_o_db);
xlabel('$B_s$ [Hz]');
ylabel('$\xi$ [dB]');
utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig3,[prefix '3']);

%% Figure 5, Cross-Correlator SNR

snr_i_db = -30:.1:30;
snr_i_lin = 10.^(snr_i_db./10);

Tp = 1e-6;
BW = [1e6,1e7,1e8,1e9]';
tbwp_db = 10*log10(Tp*BW);
snr_o_lin = Tp*BW*snr_i_lin.^2./(1+2*snr_i_lin);
snr_o_db = 10*log10(snr_o_lin);

snr_o_ideal = snr_i_db + tbwp_db;


fig5=figure;
plot(snr_i_db,snr_o_db);
hold on;set(gca,'ColorOrderIndex',1);
plot(snr_i_db,snr_o_ideal,'-.');
legendCell = arrayfun(@(x) sprintf('TB = %d ',x),Tp*BW,'UniformOutput',false);
legend(legendCell,'Location','NorthWest');
xlabel('$\xi_i$ [dB]');
ylabel('$\xi_o$ [dB]');
utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig5,[prefix '5']);

%% Figure 6 - Comparison of Performance

% Vary Time-Bandwidth Product
TBWP_db = 10:10:30;
TBWP_lin = 10.^(TBWP_db/10);

xi_in_db = -20:.1:10;
xi_in_lin = 10.^(xi_in_db/10);

[MM,XI] = meshgrid(TBWP_lin,xi_in_lin);

xi_out_lin = MM.*XI.^2./(1+2*XI);
xi_out_db = 10*log10(xi_out_lin);
XI_out = xi_out_lin;

% Energy Detector Performance
PFA = 1e-6;

if use_stat_toolbox
    eta = chi2inv(1-PFA,2*MM);
    PD = 1-ncx2cdf(eta,2*MM,2*MM.*XI);
else
    eta = utils.chi2inv(1-PFA,2*MM);
    PD = 1-utils.ncx2cdf(eta,2*MM,2*MM.*XI);
end

% legendCell = arrayfun(@(x) sprintf('ED, TB=%d dB',x),TBWP_db,'UniformOutput',false);

fig6=figure;hold on;
% for i=1:numel(legendCell)
%     set(gca,'ColorOrderIndex',i)
    h=plot(xi_in_db,PD,'DisplayName','ED');
% end
utils.excludeFromLegend(h(2:end));

% Cross-Correlator Performance
if use_stat_toolbox
    eta = chi2inv(1-PFA,2);
    PD2 = 1-ncx2cdf(eta./(1+2*XI),2,2*XI_out);
else
    eta = utils.chi2inv(1-PFA,2);
    PD2 = 1-utils.ncx2cdf(eta./(1+2*XI),2,2*XI_out);
end

hold on;
set(gca,'ColorOrderIndex',1);
h=plot(xi_in_db,PD2,'--','DisplayName','XC');
utils.excludeFromLegend(h(2:end));

% Monte Carlo Trial
xi_in_mc_db = xi_in_db(1:10:end);
xi_in_mc = 10.^(xi_in_mc_db/10);

nMC = 1e4;
nM = numel(TBWP_lin);
nXi = numel(xi_in_mc);

% Generate noise vectors
sn2 = 1; % Unit Variance

det_ed = zeros(nXi,nM);
det_xc = zeros(nXi,nM);
for ii=1:nM
    % Generate the noise vectors
    n1 = sqrt(sn2/2)*(randn(TBWP_lin(ii),nMC)+1i*randn(TBWP_lin(ii),nMC));
    n2 = sqrt(sn2/2)*(randn(TBWP_lin(ii),nMC)+1i*randn(TBWP_lin(ii),nMC));
    
    % Generate a signal vector
    s0 = sqrt(1/2)*(randn(TBWP_lin(ii),nMC)+1i*randn(TBWP_lin(ii),nMC));
    phi = exp(1i*rand(1,nMC)*2*pi);
    for jj=1:nXi
      % Scale the signal power to match SNR
      s = s0 * sqrt(xi_in_mc(jj));
      
      y1 = s+n1;
      y2 = s.*phi+n2;
      
      detResult = detector.squareLaw(y1,sn2/2,PFA);
      det_ed(jj,ii) = sum(detResult)/nMC;
          
      detResult = detector.xcorr(y1,y2,sn2,TBWP_lin(ii),PFA);
      det_xc(jj,ii) = sum(detResult)/nMC;
      
    end
end


set(gca,'ColorOrderIndex',1);
h=plot(xi_in_mc_db,det_ed,'^','DisplayName','ED (Monte Carlo)');
utils.excludeFromLegend(h(2:end));
set(gca,'ColorOrderIndex',1);
h=plot(xi_in_mc_db,det_xc,'x','DisplayName','XC (Monte Carlo)');
utils.excludeFromLegend(h(2:end));
legend('Location','SouthEast');

% Create ellipse -- 30 dB
cmap = get(gca,'ColorOrder');
annotation(fig6,'ellipse',...
    [0.377736111111111 0.594135802469136 0.0910138888888889 0.052469135802469],...
    'Color',cmap(3,:));

% Create ellipse
annotation(fig6,'ellipse',...
    [0.530513888888888 0.520061728395061 0.0910138888888889 0.052469135802469],...
    'Color',cmap(2,:));

% Create ellipse
annotation(fig6,'ellipse',...
    [0.66679861111111 0.473765432098765 0.14570138888889 0.0524691358024691],...
    'Color',cmap(1,:));

text(-15,.6,'TB=1000','FontSize',10,'Color',cmap(3,:));
hh=line([-11.8 -10.3],[.6 .58],'LineWidth',.5,'Color',cmap(3,:));
utils.excludeFromLegend(hh);
text(-15,.5,'TB=100','FontSize',10,'Color',cmap(2,:));
hh=line([-11.8 -5.4],[.5 .5],'LineWidth',.5,'Color',cmap(2,:));
utils.excludeFromLegend(hh);
text(-15,.4,'TB=10','FontSize',10,'Color',cmap(1,:));
hh=line([-11.8 -1],[.4 .44],'LineWidth',.5,'Color',cmap(1,:));
utils.excludeFromLegend(hh);
xlabel('$\xi_i$ [dB]');ylabel('$P_{\mathrm{D}}$');
grid off;

% Output the result
utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig6,[prefix '6']);

%% Figure 7 - SNR as a function of Range
fig7 = ex4_1;
utils.exportPlot(fig7,[prefix '7']);

%% Figure 8 - SNR as a function of range
fig8 = ex4_2;
utils.exportPlot(fig8,[prefix '8']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;