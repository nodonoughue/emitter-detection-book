% Draw Figures - Appendix D
%
% This script generates all of the figures that appear in
% Appendix D of the textbook.
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
prefix = fullfile(dirNm,'figD_');

% Initialize Plot Preference
utils.initPlotSettings;

% Reset the random number generator, to ensure reproducability
rng('default') ; 

%% Figure 1 - Noise vs. Noise Temp
t_ext = 0:300;
T_total = utils.constants.T0 + t_ext;

bw = 1e3;
N = noise.thermal_noise(bw(:),[],t_ext);
N_0 = noise.thermal_noise(bw(:));

fig1=figure;
semilogx(t_ext,N-N_0);
xlabel('Combined Sky Noise Temperature [K]');
ylabel('Increase in Noise Level [dB]');
grid on;

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig1,[prefix '1']);

%% Figure 2 - Cosmic Noise

% Plot cosmic noise [dB] as a function of frequency for a fixed bandwidth
% (a) without solar/lunar gain
% (b) with solar gain = 0 dBi
% (c) with solar gain = 30 dBi
% (d) with lunar gain = 0 dBi
% (e) with lunar gain = 30 dBi

freq = 100e6:100e6:1e10;
T_a = noise.cosmic_noise_temp(freq,0,.95);
T_b = noise.cosmic_noise_temp(freq,0,.95,[],0);
T_c = noise.cosmic_noise_temp(freq,0,.95,[],30);
T_d = noise.cosmic_noise_temp(freq,0,.95,0);
T_e = noise.cosmic_noise_temp(freq,0,.95,30);

fig2=figure;
loglog(freq/1e9,T_a,'-','LineWidth',1,'DisplayName','Cosmic Noise');
hold on;
set(gca,'ColorOrderIndex',3);
hh=loglog(freq/1e9,T_c,'LineWidth',1,'DisplayName','Cosmic with G_m=30 dBi');
utils.excludeFromLegend(hh);
set(gca,'ColorOrderIndex',4);
hh=loglog(freq/1e9,T_e,'LineWidth',1,'DisplayName','Cosmic with G_s=30 dBi');
utils.excludeFromLegend(hh);
loglog(freq/1e9,utils.constants.T0*ones(size(freq)),'k:','LineWidth',1,'DisplayName','Thermal Noise');
set(gca,'ColorOrderIndex',1);
loglog(freq/1e9,T_a+utils.constants.T0,'-.','LineWidth',1,'DisplayName','Thermal + Cosmic Noise');
set(gca,'ColorOrderIndex',3);
hh=loglog(freq/1e9,T_c+utils.constants.T0,'-.','LineWidth',1,'DisplayName','Cosmic with G_m=30 dBi');
utils.excludeFromLegend(hh);
set(gca,'ColorOrderIndex',4);
hh=loglog(freq/1e9,T_e+utils.constants.T0,'-.','LineWidth',1,'DisplayName','Cosmic with G_s=30 dBi');
utils.excludeFromLegend(hh);

annotation(fig2,'arrow',[0.582465277777775 0.582465277777775],...
    [0.658489583333333 0.801302083333333]);
text(1,200,'Impact of cosmic noise','FontSize',9);
text(.65,3,'Sidelobe Cosmic Noise','FontSize',9);
text(1.5,9,'Mainbeam pointed at Moon','FontSize',9);
text(2.4,75,'Mainbeam pointed at Sun','FontSize',9);
legend('Location','NorthWest');
xlabel('Frequency [GHz]');
ylabel('Noise Temperature [K]');

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig2,[prefix '2']);

%% Figure 3 - Atmospheric Noise
zenith_angle_deg = [0,10,30,60];
zenith_angle_rad = zenith_angle_deg*pi/180;

% Set up frequencies
[fo,fw] = atm.getSpectralLines();
f = sort([fo,fw,fo+50e6,fw+50e6,fo-100e6,fw-100e6,1e9:1e9:350e9],'ascend');

Ta = zeros(numel(f),numel(zenith_angle_deg));
for idx_ang = 1:numel(zenith_angle_deg)
    Ta(:,idx_ang) = noise.atmospheric_noise(f(:),0,90-zenith_angle_deg(idx_ang));
end

fig3=figure;
semilogx(f/1e9,Ta);
hold on;
semilogx(f/1e9,utils.constants.T0*ones(size(f)),'k:','DisplayName','Thermal Noise');
set(gca,'ColorOrderIndex',1);
hh=semilogx(f/1e9,Ta+utils.constants.T0,'-.','DisplayName','Thermal + Atmospheric Noise');
utils.excludeFromLegend(hh(2:end));
legend(cat(2,arrayfun(@(x) sprintf('%d%s from Zenith',x,char(176)),zenith_angle_deg,'UniformOutput',false),'Noise Temperature','Atmospheric + Noise'),'Location','NorthWest');
xlabel('Freq [GHz]');
ylabel('Noise Temperature[K]');
xlim([1 350]);

annotation(fig3,'arrow',[0.72 0.72],...
    [0.547839506172839 0.861111111111111]);
text(60,310,'Impact of Atmospheric Noise','FontSize',9);
utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig3,[prefix '3']);

%% Figure 4 - Ground Noise
Gg = -30:0;

T = noise.ground_noise(Gg);

fig4=figure;
plot(Gg,T,'DisplayName','Ground Noise');
hold on;
plot(Gg,utils.constants.T0*ones(size(Gg)),'k:','DisplayName','Thermal Noise');
plot(Gg,T+utils.constants.T0,'k-.','DisplayName','Thermal + Ground Noise');
xlabel('Average Ground Antenna Gain [dBi]');
ylabel('Noise Temperature [K]');
legend('Location','NorthWest');
% set(gca,'yscale','log');

annotation(fig4,'arrow',[0.91 0.91],...
    [0.74437037037037 0.819444444444444]);
text(-8,270,'Impact of Ground Noise','FontSize',9);
utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig4,[prefix '4']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;