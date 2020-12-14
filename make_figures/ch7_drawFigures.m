% Draw Figures - Chapter 7
%
% This script generates all of the figures that appear in
% Chapter 7 of the textbook.
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

prefix = fullfile(dirNm,'fig7_');

% Initialize Plot Preference
utils.initPlotSettings;

% Initializes colorSet - Mx3 RGB vector for successive plot lines
colors = get(groot,'DefaultAxesColorOrder');

% Reset the random number generator, to ensure reproducability
rng('default') ; 

% Add path to folder for examples from the textbook
addpath('examples');

%% Figure 1 - Lobing Interpolation
numLobes = 3;
numSamples = 501;
x = linspace(-numLobes,numLobes,numSamples);

x0 = .3;

y = sinc(x-x0);

% Sample points
xs = -2:.5:2;
ys = sinc(xs-x0);

th_mult = 45/numLobes; % Set the end points to +/- 45 degrees

fig1=figure;
plot(xs*th_mult,abs(ys),'ko','DisplayName','Sample Points');
hold on;
plot(x*th_mult,abs(y),'k--','DisplayName','Interpolated Beampattern');
plot([1 1]*x0*th_mult,[0 1],'k:','DisplayName','True Emitter Bearing');
xlabel('$\phi$ [degrees]');
ylabel('Antenna Gain');
legend('Location','NorthEast');

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig1,[prefix '1']);

%% Figure 3 -- Adcock CRLB 
if force_recalc
% Create the antenna pattern generating function
d_lam = .25;
[g,g_dot] = aoa.make_gain_functions('adcock',d_lam,0);
    % --- NOTE --- g,g_dot take radian inputs (psi, not theta)
    
% Generate the angular samples and true gain values
th_true = 5;
psi_true = th_true*pi/180;
psi_res = .001; % desired resolution from multi-stage search;
                % see aoa.directional_df for details.
N = 10;
th = linspace(-180,180-360/N,N); % evenly spaced across -45:45 degrees
psi = th(:)*pi/180;
x = g((psi-psi_true)); % Actual gain values

% Set up the parameter sweep
M_vec = [1,10,100]; % Number of temporal samples at each antenna test point
snr_db_vec = -20:2:20; % signal to noise ratio
num_mc = 1000; % number of monte carlo trials at each parameter setting

% Set up output scripts
rmse_psi = zeros(numel(M_vec),numel(snr_db_vec));
crlb_psi = zeros(numel(M_vec),numel(snr_db_vec));

% Loop over parameters
fprintf('Executing Adcock Monte Carlo sweep...\r\n');
for idx_M = 1:numel(M_vec)
    M = M_vec(idx_M);
    this_num_mc = num_mc;
    fprintf('\t M=%d',M);

    % Generate Monte Carlo Noise with unit power
    % -- for simplicity, we only generate the real component, since this
    %    receiver is only working on the real portion of the received
    %    signal
    noise_base = (1/sqrt(2))*randn(N,M,this_num_mc);
    for idx_snr = 1:numel(snr_db_vec)
        fprintf('.');
        
        % Compute noise power, scale base noise
        snr_db = snr_db_vec(idx_snr);
        noise_pwr = 10.^(-snr_db/10);
        
        % Generate noisy measurement
        n = sqrt(noise_pwr)*noise_base;
        y = x+n;
        
        % Estimate
        psi_est = zeros(1,this_num_mc);
        for idx_mc = 1:this_num_mc
            psi_est(idx_mc) = aoa.directional_df(y(:,:,idx_mc),psi,g,psi_res,min(psi),max(psi));
        end
        rmse_psi(idx_M,idx_snr) = sqrt(sum(abs((psi_est-psi_true)).^2)/this_num_mc);
        
        % CRLB
        crlb_psi(idx_M,idx_snr) = abs(aoa.directional_crlb(snr_db,M,g,g_dot,psi,psi_true));
    end
fprintf('done.\n');
end


fig3=figure;
for idx_M=1:numel(M_vec)
    set(gca,'ColorOrderIndex',idx_M);
    h1=semilogy(snr_db_vec,sqrt(crlb_psi(idx_M,:))*180/pi,'DisplayName','CRLB');%sprintf('CRLB, M=%d',M_vec(idx_M)));
    hold on;
    set(gca,'ColorOrderIndex',idx_M);
    h2=semilogy(snr_db_vec,rmse_psi(idx_M,:)*180/pi,'--','DisplayName','Simulation Result');%sprintf('Simulation Result, M=%d',M_vec(idx_M)));
    if idx_M~=1
        utils.excludeFromLegend(h1);
        utils.excludeFromLegend(h2);
    end
end
xlabel('$\xi$ [dB]');
ylabel('RMSE [deg]');
title('Adcock DF Performance');
legend('Location','SouthWest');

text(2,17,'M=1','FontSize',10);
text(2,4,'M=10','FontSize',10);
text(2,1.4,'M=100','FontSize',10);

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig3,[prefix '3']);
end


%% Figure 5 -- Rectangular Aperture CRLB example
if force_recalc
% Create the antenna pattern generating function
D_lam = 5;
[g,g_dot] = aoa.make_gain_functions('rectangular',D_lam,0);
    % --- NOTE --- g,g_dot take radian inputs (psi, not theta)
    
% Generate the angular samples and true gain values
th_true = 5;
psi_true = th_true*pi/180;
N = 36; % number of samples
th = linspace(-180,180-360/N,N);
psi = th(:)*pi/180;
psi_res = .001; % desired resolution from multi-stage search;
                % see aoa.directional_df for details.
x = g((psi-psi_true)); % Actual gain values

% Set up the parameter sweep
M_vec = [1,10,100]; % Number of temporal samples at each antenna test point
snr_db_vec = -10:20; % signal to noise ratio
num_mc = 10000; % number of monte carlo trials at each parameter setting

% Set up output scripts
rmse_psi = zeros(numel(M_vec),numel(snr_db_vec));
crlb_psi = zeros(numel(M_vec),numel(snr_db_vec));

% Loop over parameters
fprintf('Executing Rectangular Aperture Monte Carlo sweep...\r\n');
for idx_M = 1:numel(M_vec)
    M = M_vec(idx_M);
    this_num_mc = num_mc / M;
    fprintf('\t M=%d',M);
    
    % Generate Monte Carlo Noise with unit power
    noise_base = 1/sqrt(2)*(randn(N,M,this_num_mc)+1i*randn(N,M,this_num_mc));

    for idx_snr = 1:numel(snr_db_vec)
        fprintf('.');
        
        % Compute noise power, scale base noise
        snr_db = snr_db_vec(idx_snr);
        noise_pwr = 10.^(-snr_db/10);
        
        % Generate noisy measurement
        n = sqrt(noise_pwr)*noise_base;
        y = x+n;
        
        % Estimate
        psi_est = zeros(1,this_num_mc);
        for idx_mc = 1:this_num_mc
            psi_est(idx_mc) = aoa.directional_df(y(:,:,idx_mc),psi,g,psi_res,-pi,pi);
        end
        rmse_psi(idx_M,idx_snr) = sqrt(sum(abs((psi_est-psi_true)).^2)/this_num_mc);
        
        % CRLB
        crlb_psi(idx_M,idx_snr) = abs(aoa.directional_crlb(snr_db,M,g,g_dot,psi,psi_true));
    end
fprintf('done.\n');
end


fig5=figure;
for idx_M=1:numel(M_vec)
    set(gca,'ColorOrderIndex',idx_M);
    h1=semilogy(snr_db_vec,sqrt(crlb_psi(idx_M,:))*180/pi,'DisplayName','CRLB');%sprintf('CRLB, M=%d',M_vec(idx_M)));
    hold on;
    set(gca,'ColorOrderIndex',idx_M);
    h2=semilogy(snr_db_vec,rmse_psi(idx_M,:)*180/pi,'--','DisplayName','Simulation Result');%sprintf('Simulation Result, M=%d',M_vec(idx_M)));
    if idx_M~=1
        utils.excludeFromLegend(h1);
        utils.excludeFromLegend(h2);
    end
end
xlabel('$\xi$ [dB]');
ylabel('RMSE [deg]');
title('Rectangular Aperture DF Performance');
legend('Location','SouthWest');

text(7,5,'M=1','FontSize',10);
text(7,1.7,'M=10','FontSize',10);
text(7,.5,'M=100','FontSize',10);

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig5,[prefix '5']);

end

%% Figure 6b -- Watson-Watt Patterns

phi = -pi:pi/101:pi;

G1 = 2*cos(phi);
G2 = 2*cos(phi-pi/2);
Gomni = ones(size(phi));

fig6b=figure;
p1=polar(phi,abs(G1));hold on;
p2=polar(phi,abs(G2));
p3=polar(phi,Gomni,'--');

% find all of the lines in the polar plot
h = findall(gcf,'type','line');
% remove the handle for the polar plot line from the array
h(h == p1) = [];
h(h == p2) = [];
h(h == p3) = [];
% delete all other lines
delete(h);
h = findall(gcf,'type','Patch');
delete(h);
% If you want, use the following lines to remove the text.
% (I prefer to leave it in place)
% find and remove all of the text objects in the polar plot
delete(findall(gcf,'type','text'))
text(1.1,.25,'Reference','FontSize',10);
text(1.5,1,'Horizontal Adcock','FontSize',10);
text(0,2.1,'Vertical Adcock','FontSize',10);


utils.exportPlot(fig6b,[prefix,'6b']);

%% Figure 7 -- Watson Watt Performance
if force_recalc
    
% Generate the Signals
th_true = 45;
psi_true = th_true*pi/180;
f = 1e9;
t_samp = 1/(3*f); % ensure the Nyquist criteria is satisfied    

% Set up the parameter sweep
M_vec = [1,10,100]; % Number of temporal samples at each antenna test point
snr_db_vec = -10:.2:20; % signal to noise ratio
num_mc = 1e5; % number of monte carlo trials at each parameter setting

% Set up output scripts
rmse_psi = zeros(numel(M_vec),numel(snr_db_vec));
crlb_psi = zeros(numel(M_vec),numel(snr_db_vec));

% Loop over parameters
fprintf('Executing Watson-Watt Monte Carlo sweep...\r\n');
for idx_M = 1:numel(M_vec)
    M = M_vec(idx_M);
    this_num_mc = num_mc/M;
    fprintf('\t M=%d',M);
    
    % Generate signal vectors
    t_vec = (0:M-1)*t_samp;
    r0 = cos(2*pi*f*t_vec).';
    y0 = sind(th_true)*r0;
    x0 = cosd(th_true)*r0;
    
    % Generate Monte Carlo Noise with unit power
    noise_base_r = sqrt(sum(abs(r0).^2)/M)*randn(M,this_num_mc);
    noise_base_x = sqrt(sum(abs(r0).^2)/M)*randn(M,this_num_mc);
    noise_base_y = sqrt(sum(abs(r0).^2)/M)*randn(M,this_num_mc);

    for idx_snr = 1:numel(snr_db_vec)
        if mod(idx_snr,10)==0
            fprintf('.');
        end
        
        % Compute noise power, scale base noise
        snr_db = snr_db_vec(idx_snr);
        noise_pwr = 10.^(-snr_db/10);
        
        % Generate noisy measurement
        r = r0 + noise_base_r*sqrt(noise_pwr);
        y = y0 + noise_base_y*sqrt(noise_pwr);
        x = x0 + noise_base_x*sqrt(noise_pwr);
        
        % Compute the estimate
        psi_est = zeros(1,this_num_mc);
        for idx_mc = 1:this_num_mc
            psi_est(idx_mc) = aoa.watson_watt_df(r(:,idx_mc),x(:,idx_mc),y(:,idx_mc));
        end

        rmse_psi(idx_M,idx_snr) = sqrt(sum(abs((psi_est-psi_true)).^2)/this_num_mc);
        
        % CRLB
        crlb_psi(idx_M,idx_snr) = abs(aoa.watson_watt_crlb(snr_db,M));
    end
fprintf('done.\n');
end


fig7=figure;
for idx_M=1:numel(M_vec)
    set(gca,'ColorOrderIndex',idx_M);
    h1=semilogy(snr_db_vec,sqrt(crlb_psi(idx_M,:))*180/pi,'DisplayName','CRLB');%sprintf('CRLB, M=%d',M_vec(idx_M)));
    hold on;
    set(gca,'ColorOrderIndex',idx_M);
    h2=semilogy(snr_db_vec,rmse_psi(idx_M,:)*180/pi,'--','DisplayName','Simulation Result');%sprintf('Simulation Result, M=%d',M_vec(idx_M)));
    if idx_M~=1
        utils.excludeFromLegend(h1);
        utils.excludeFromLegend(h2);
    end
end
xlabel('$\xi$ [dB]');
ylabel('RMSE [deg]');
title('Watson-Watt DF Performance');
legend('Location','SouthWest');

text(10,25,'M=1','FontSize',10);
text(10,7,'M=10','FontSize',10);
text(10,2.5,'M=100','FontSize',10);

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig7,[prefix '7']);

end

%% Figure 8b -- Doppler
M=128;
phi0=5*pi/3;
fd = cos((1:M)/M*2*pi+phi0);
f0 = 2.5;

fig8b=figure;
plot(1:M,f0*ones(size(fd)),'DisplayName','f');hold on;
plot(1:M,f0+fd,'DisplayName','f_d');
xlabel('Time');
ylabel('Frequency');
% legend('Location','NorthWest');
text(.8*M,f0+.1,'$f$','FontSize',12);
text(.35*M,f0+fd(fix(.3*M))-.15,'$f+f_d(t)$','FontSize',12);
fd_desc = fd;
fd_desc(fd(2:end)>fd(1:end-1))=Inf;
[~,idx] = min(abs(fd_desc));
text(idx+1,f0-.8*max(fd),'$\tau$','FontSize',12);
plot(idx*[1 1],f0+[min(fd) 0],'k--');
utils.setPlotStyle(gca,{'box only','widescreen'});
utils.exportPlot(fig8b,[prefix,'8b']);

%% Figure 10 -- Doppler CRLB
if force_recalc
    
% Generate the Signals
th_true = 45;
psi_true = th_true*pi/180;
A = 1;
phi0 = 2*pi*rand();
f = 1e9;
ts = 1/(5*f);

% Doppler antenna parameters
psi_0 = 0; % initial pointing angle [radians]
c=3e8;
lam = c/f;
R = lam/2;
psi_res = .0001; % desired doppler resolution

% Set up the parameter sweep
M_vec = [10,100,1000]; % Number of temporal samples at each antenna test point
snr_db_vec = -10:2:30; % signal to noise ratio
num_mc = 100000; % number of monte carlo trials at each parameter setting

% Set up output scripts
rmse_psi = zeros(numel(M_vec),numel(snr_db_vec));
crlb_psi = zeros(numel(M_vec),numel(snr_db_vec));

% Loop over parameters
fprintf('Executing Doppler Monte Carlo sweep...\r\n');
for idx_M = 1:numel(M_vec)
    M = M_vec(idx_M);
    this_num_mc = num_mc;
    fprintf('\t M=%d',M);
    
    % Reference signal
    t_vec = ts*(0:M-1);
    r0 = A*exp(1i*phi0)*exp(1i*2*pi*f*t_vec);
    
    % Doppler signal
    fr = 1/(ts*M); % Ensure a single cycle during M
    x0 = A*exp(1i*phi0)*exp(1i*2*pi*f*t_vec).*exp(1i*2*pi*f*R/c*cos(2*pi*fr*t_vec-psi_true));
    
    % Generate noise signal
    noise_base_r = sqrt(A/2)*(randn(M,this_num_mc)+1i*randn(M,this_num_mc));
    noise_base_x = sqrt(A/2)*(randn(M,this_num_mc)+1i*randn(M,this_num_mc));

    for idx_snr = 1:numel(snr_db_vec)
%         if mod(idx_snr,10)==0
            fprintf('.');
%         end
        
        % Compute noise power, scale base noise
        snr_db = snr_db_vec(idx_snr);
        noise_pwr = 10.^(-snr_db/10);
        
        % Generate noisy signals
        r = r0(:) + noise_base_r*sqrt(noise_pwr);
        x = x0(:) + noise_base_x*sqrt(noise_pwr);
        
        % Compute the estimate
        psi_est = zeros(1,this_num_mc);
        for idx_mc = 1:this_num_mc
            psi_est(idx_mc) = aoa.doppler_df(r(:,idx_mc),x(:,idx_mc),ts,f,R,fr,psi_res,-pi,pi);
        end
        rmse_psi(idx_M,idx_snr) = sqrt(sum(abs((psi_est-psi_true)).^2)/this_num_mc);
        
        % CRLB
        crlb_psi(idx_M,idx_snr) = aoa.doppler_crlb(snr_db,M,A,ts,f,R,fr,psi_true);
    end
fprintf('done.\n');
end


fig10=figure;
for idx_M=1:numel(M_vec)
    set(gca,'ColorOrderIndex',idx_M);
    h1=semilogy(snr_db_vec,sqrt(crlb_psi(idx_M,:))*180/pi,'DisplayName','CRLB');%sprintf('CRLB, M=%d',M_vec(idx_M)));
    hold on;
    set(gca,'ColorOrderIndex',idx_M);
    h2=semilogy(snr_db_vec,rmse_psi(idx_M,:)*180/pi,'--','DisplayName','Simulation Result');%sprintf('Simulation Result, M=%d',M_vec(idx_M)));
    if idx_M~=1
        utils.excludeFromLegend(h1);
        utils.excludeFromLegend(h2);
    end
end
xlabel('$\xi$ [dB]');
ylabel('RMSE [deg]');
title('Doppler DF Performance');
legend('Location','SouthWest');

text(5,5,'M=10','FontSize',10);
text(5,1.7,'M=100','FontSize',10);
text(5,.5,'M=1,000','FontSize',10);

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig10,[prefix '10']);

end

%% Figure 12 -- Interferometer

d_lam = [.5,1,2,4];

psi_s = linspace(-60,60,201)*pi/180;
psi_0 = -10*pi/180;

g_true = (1+exp(1i*2*pi*d_lam(:)*(sin(psi_s)-sin(psi_0))))/2;

fig12=figure;hold on;
for idx_plot = 1:numel(d_lam)
    this_width = 2 - idx_plot*.25;
    plot(180*psi_s/pi,10*log10(abs(g_true(idx_plot,:))),'linewidth',this_width);
end
%plot(180*th_s/pi,10*log10(abs(g_true)));
ylim([-10 0]);
xlabel('$\theta$');
ylabel('Normalized Response [dB]');
legend(arrayfun(@(x) sprintf('d/\\lambda=%.1f',x),d_lam,'UniformOutput',false));
hold on;
plot(180*psi_0/pi*[1 1],[-10 0],'k--','DisplayName','AOA');

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig12,[prefix '11']);

%% Figure 14 -- Interferometer Example
fig14 = ex7_1;
ex7_2(fig14);
ex7_3(fig14);
ex7_4(fig14);

utils.setPlotStyle(gca,{'widescreen','tight'});
%utils.exportPlot(fig14,[prefix '14']);

%% Figure 15b -- Monopulse
% Make monopulse beampatterns

N = 10001;
x = linspace(-3,3,N);

% Illumination Pattern
a = sinc(x);

% Individual Signals
s = a.*cos(pi*x);
d = a.*sin(pi*x);
k = d./s;
k(abs(x)>.5)=Inf;

fig15b=figure;
plot(x,s,'DisplayName','Sum');
hold on;
set(gca,'ColorOrderIndex',3);
plot(x,d,'DisplayName','Difference');
plot(x,k,'--','DisplayName','Monopulse Ratio');
xlabel('Angle (beamwidths)');
ylabel('Normalized Antenna Pattern');
legend('Location','SouthWest');
ylim([-1 1]);

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig15b,[prefix '15b']);

%% Figure 16 -- Monopulse Error

th_bw = 4;
th = -th_bw/2:.1:th_bw/2;
km = 1.606;

snr_db = [10 20]';
snr_lin = 10.^(snr_db/10);

rmse_th = sqrt(1./(km^2.*snr_lin).*(1+(km*th/th_bw).^2));

fig16=figure;hold on;
plot(th/th_bw,rmse_th(1,:));
set(gca,'ColorOrderIndex',3);
plot(th/th_bw,rmse_th(2,:));
legend('\xi = 10 dB','\xi = 20 dB');
xlabel('$\overline{\theta}$');
ylabel('$\sigma_{\overline\theta}$ [deg]');
ylim([0 .5]);
grid on;

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig16,[prefix '16']);

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;