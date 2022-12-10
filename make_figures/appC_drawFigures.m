% Draw Figures - Appendix C
%
% This script generates all of the figures that appear in
% Appendix C of the textbook.
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
prefix = fullfile(dirNm,'figC_');

% Initialize Plot Preference
utils.initPlotSettings;

% Reset the random number generator, to ensure reproducability
rng('default') ; 

%% Figure 2 - Dry Air and Water Vapor
alt=[0,10,20]*1e3;
atmStruct = atm.standardAtmosphere(alt);

% Set up frequencies
[fo,fw] = atm.getSpectralLines();
f = sort([fo,fw,fo+50e6,fw+50e6,fo-100e6,fw-100e6,1e9:1e9:350e9],'ascend');

% Compute Loss Coeffs
[ao,aw] = atm.gasLossCoeff(f(:),atmStruct.P,atmStruct.e,atmStruct.T);

fig2=figure;
h=loglog(f/1e9,ao,':','DisplayName','Dry Air Only');hold on;
utils.excludeFromLegend(h(2:end));
set(gca,'ColorOrderIndex',1);
h=loglog(f/1e9,aw,'--','DisplayName','Water Vapor');
utils.excludeFromLegend(h(2:end));
set(gca,'ColorOrderIndex',1);
h=loglog(f/1e9,ao+aw,'-','DisplayName','Total');
utils.excludeFromLegend(h(2:end));
xlim([1 350]);
xlabel('Frequency [GHz]');
ylabel('Gas Loss Coefficient $\gamma_g$ [dB/km]');
legend;%('\gamma_o','\gamma_w','Total');
ylim([1e-5,1e2]);
text(2,1.5e-2,'0 km');
text(2,1.5e-3,'10 km');
text(2,8e-5,'20 km');
legend('Location','NorthWest');
utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig2,[prefix '2']);

%% Figure 3 - Rain Loss Coefficient
rainRateSet = [1,4,16,100];
rainRateNames = {'Light Rain','Moderate Rain','Heavy Rain','Very Heavy Rain'};

freq_ghz = 1:.5:100;
el_ang = 0*pi/180;
pol_ang = [0 pi/2];
pol_set = {'Horizontal','Vertical'};

gamma = atm.rainLossCoeff(reshape(freq_ghz,1,1,[])*1e9,pol_ang(:),el_ang,rainRateSet(:)');

fig3=figure;
hold on;
set(gca,'LineStyleOrder',{'-','-.'});
for ii=1:numel(pol_ang)
    set(gca,'LineStyleOrderIndex',ii);
    for jj=1:numel(rainRateSet)
        set(gca,'ColorOrderIndex',jj);
        %thisName = rainRateNames{jj};
        thisName = pol_set{ii};
        
        hh=plot(freq_ghz,squeeze(gamma(ii,jj,:)),'DisplayName',thisName);
        if jj>1
            utils.excludeFromLegend(hh);
        end
    end
end
set(gca,'yscale','log');
set(gca,'xscale','log');
grid on
ylim([.01 50])
% Add rainfall condition labels
ht=text(10,6,'Very Heavy');set(ht,'rotation',25);
ht=text(10,.6,'Heavy');set(ht,'rotation',40);
ht=text(10,.12,'Moderate');set(ht,'rotation',40);
ht=text(10,.023,'Light');set(ht,'rotation',40);

xlabel('Frequency [GHz]');
ylabel('Rain Loss Coefficient $\gamma_r$ [dB/km]');
legend('Location','NorthWest');

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig3,[prefix '3']);


%% Figure 4
fogSet = [.032, .32, 2.3];
fogNames = {'600 m Visibility','120 m Visibility','30 m Visibility'};

freq_ghz = 1:.5:100;

gamma = atm.fogLossCoeff(freq_ghz*1e9,fogSet(:));

fig4=figure;
hold on;
for ii=1:numel(fogSet)
    set(gca,'ColorOrderIndex',1+2*(ii-1));
    plot(freq_ghz,gamma(ii,:),'DisplayName',fogNames{ii});
end
set(gca,'xscale','log');
set(gca,'yscale','log');
grid on;
ht=text(10,.2,'30 m Visibility');set(ht,'rotation',30);
ht=text(10,.025,'120 m Visibility');set(ht,'rotation',30);
ht=text(32,.025,'600 m Visibility');set(ht,'rotation',30);

ylim([.01 10]);
xlabel('Frequency [GHz]');
ylabel('Cloud Loss Coefficient $\gamma_c$ [dB/km]');
% legend('Location','NorthWest');

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig4,[prefix '4']);


%% Figure 5 -- Zenith Loss

% To compute zenith loss, we step through different altitudes, and
% accumulate the loss at each altitude
alt_step = 100;
max_alt = 100e3;
alt_lower_vec = 0:alt_step:(max_alt-alt_step);
num_alt = numel(alt_lower_vec);

% Frequencies over which to calculate
% Set up frequencies
[fo,fw] = atm.getSpectralLines();
f = sort([fo,fw,fo+100e6,fw+100e6,fo-100e6,fw-100e6,1e9:1e9:350e9],'ascend');

nadir_deg_set = [0,10,30,60];
legend_entries = arrayfun(@(x) sprintf('%d%s from Zenith',x,char(176)), nadir_deg_set,'UniformOutput',false);
legend_entries{1} = 'Zenith';
loss = atm.calcZenithLoss(f(:),0,nadir_deg_set);

fig5=figure;
set(gca,'LineStyleOrder','-|--|-.|:');
set(gca,'ColorOrder',[0 0 0]);
set(gca,'yscale','log');
set(gca,'xscale','log');
hold on;
loglog(f/1e9,loss,'DisplayName','Zenith');
legend(legend_entries);

xlabel('Frequency [GHz]');
legend('Location','NorthWest');
ylabel('Zenith Attenuation [linear]');
grid on;
xlim([1 350])

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig5,[prefix '5']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;