% Draw Figures - Appendix B
%
% This script generates all of the figures that appear in
% Appendix B of the textbook.
%
% Nicholas O'Donoughue
% 1 July 2019

% Flag to force re-execution of long scripts
if ~exist('force_recalc','var')
    force_recalc = false;
end

% Clear Figures
close all;

% Set up directory and filename for figures
dirNm = fullfile(pwd,'figures');
dirNmComponents = fullfile(pwd,'Graphics','Components');
if ~exist(dirNm,'dir')
    mkdir(dirNm);
end
if ~exist(dirNmComponents,'dir')
    mkdir(dirNmComponents);
end

prefix = fullfile(dirNm,'figB_');

% Initialize Plot Preference
utils.initPlotSettings;

% Initializes colorSet - Mx3 RGB vector for successive plot lines
colors = get(groot,'DefaultAxesColorOrder');

% Reset the random number generator, to ensure reproducability
rng('default') ; 

%% Figure 4 - Fresnel Zone Illustration

rng_m = [1e3:1e3:100e3,200e3:100e3:10000e3];

% Three Situations
% 1 - L Band, ht=hr=10 m
% 2 - L Band, ht=hr=100 m
% 3 - X Band, ht=hr=100 m

freq = [1e9,1e10];
ht = [10, 100];
hr = [10, 100];
fspl_legend = arrayfun(@(x) sprintf('Free-Space, f=%d GHz',x/1e9),freq,'UniformOutput',false);
tray_legend = arrayfun(@(x) sprintf('Two-Ray, h_T=h_R=%d m',x),ht,'UniformOutput',false);

% Path Loss
fspl = prop.freeSpacePathLoss(rng_m,freq(:),false);
tworay = prop.twoRayPathLoss(rng_m,freq(:),ht(:),hr(:),false);


fig4=figure;
set(gca,'ColorOrder',[0 0 0; .5 .5 .5]);
hold on;
plot(rng_m/1e3,fspl);
set(gca,'ColorOrderIndex',1);
plot(rng_m/1e3,tworay,'-.');
set(gca,'xscale','log')
legend(cat(2,fspl_legend,tray_legend),'Location','NorthWest');

% Overlay R_FZ

r_fz = prop.fresnelZone(freq(:),ht,hr);
y_fz = zeros(size(r_fz));
for idx_freq = 1:numel(freq)
    for idx_ht = 1:numel(ht)
        this_fz = r_fz(idx_freq,idx_ht);
        
        [~,idx_fz] = min(abs(rng_m-this_fz));
        y_fz(idx_freq,idx_ht) = fspl(idx_freq,idx_fz);
    end
end
plot(r_fz(:)/1e3,y_fz(:),'k^','DisplayName','R_{FZ}','MarkerSize',8);

xlabel('Path Length [km]');
ylabel('Loss');

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig4,[prefix '4']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;