% Draw Figures - Chapter 8
%
% This script generates all of the figures that appear in
% Chapter 8 of the textbook.
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

prefix = fullfile(dirNm,'fig8_');

% Initialize Plot Preference
utils.initPlotSettings;

% Initializes colorSet - Mx3 RGB vector for successive plot lines
colors = get(groot,'DefaultAxesColorOrder');

% Reset the random number generator, to ensure reproducability
rng('default') ; 

% Add path to folder for examples from the textbook
addpath('examples');

%% Figure 3 -- Array Factor

d_lam = 1/2;
th_0 = -30;
th = linspace(-90,90,1001);

psi_0 = th_0*pi/180;
psi = th*pi/180;

u = sin(psi);

d_vec = [10,25,100];
af = zeros(numel(th),numel(d_vec));
for idx_d = 1:numel(d_vec)
    af(:,idx_d) = array.compute_array_factor_ula(d_lam,d_vec(idx_d),psi,psi_0);
end

fig3=figure;
set(gca,'ColorOrder',[0 0 0]);hold on;
plot(u,af);
legend(arrayfun(@(x) sprintf('N= %d',x),d_vec,'UniformOutput',false));
xlabel('u');
ylabel('Array Factor [linear]');
grid off;

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig3,[prefix '3']);

%% Figure 4 -- Array Factor w/Grating Lobes
N = 11;
th_0 = -30;
th = linspace(-90,90,1001);

psi_0 = th_0*pi/180;
psi = th*pi/180;

u = sin(psi);

d_vec = [.5,1,2];
af = zeros(numel(th),numel(d_vec));
for idx_d = 1:numel(d_vec)
    af(:,idx_d) = array.compute_array_factor_ula(d_vec(idx_d),N,psi,psi_0);
end

fig4=figure;
set(gca,'ColorOrder',[0 0 0]);hold on;
plot(u,af);
legend(arrayfun(@(x) sprintf('d/\\lambda= %.1f',x),d_vec,'UniformOutput',false));
xlabel('u');
ylabel('Array Factor [linear]');
grid off;

% Annotation
hold on;
text(-.4,.8,'Mainlobe','FontSize',10);
text(.1,.65,'Grating Lobes','FontSize',10);
h1=plot([.05 .1],[.5 .6],'k-');
h2=plot([.425 .1],[.5 .6],'k-');
%h3=plot([132 139],[.8 .5],'k-');
utils.excludeFromLegend([h1 h2]);


utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig4,[prefix '4']);

%% Figure 5 - Grating Lobes w/Cosine Pattern
N = 11;
th_0 = -30;
th = linspace(-90,90,1001);

psi_0 = th_0*pi/180;
psi = th*pi/180;

u = sin(psi);

el_pat = @(psi) cos(psi).^1.2;

d_vec = [.5,1,2];
af = zeros(numel(th),numel(d_vec));
for idx_d = 1:numel(d_vec)
    af(:,idx_d) = array.compute_array_factor_ula(d_vec(idx_d),N,psi,psi_0,el_pat);
end

fig5=figure;
set(gca,'ColorOrder',[0 0 0]);
hold on;
plot(u,af);
legend(arrayfun(@(x) sprintf('d/\\lambda= %.1f',x),d_vec,'UniformOutput',false));
xlabel('u');
ylabel('Array Factor [linear]');
grid off;

% Annotation
hold on;
text(-.4,.8,'Mainlobe','FontSize',10);
text(.1,.65,'Grating Lobes','FontSize',10);
h1=plot([.05 .1],[.5 .6],'k-');
h2=plot([.425 .1],[.5 .6],'k-');
%h3=plot([132 139],[.8 .5],'k-');
utils.excludeFromLegend([h1 h2]);


utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig5,[prefix '5']);

%% Figure 6 -- Beamwidth Plot

% Set up array
N = 25;
d_lam = .5;
d_psi = .89/((N-1)*d_lam);
du = sin(d_psi);
v = array.make_steering_vector(d_lam,N);

% Set up source signals
spacing = 1;%[.8 1 1.2];
u1 = du*spacing/2;
u2 = -u1;

psi1 = asin(u1);
psi2 = asin(u2);

v1 = v(psi1); % N x numel(u1)
v2 = v(psi2); % N x numel(u2)

% Add noise
snr_db = 50;
snr = 10.^(snr_db/10);
%n = sqrt(1./(snr*2))*(randn(size(v1))+1i*randn(size(v1)));
%x = v1 + v2 + n;

% Compute Array Factor by calling with the noisy input signal
% x in the place of the beamformer h, to compute the response of a standard
% beamformer to the actual received signal.
u_scan = linspace(-.5,.5,1001);
psi_scan = asin(u_scan);
AF1 = zeros(numel(psi_scan),numel(v1));
AF2 = zeros(numel(psi_scan),numel(v1));
%AF = zeros(numel(psi_scan),numel(v1));
for i=1:numel(u1)
    AF1(:,i) = array.compute_array_factor(v,v1(:,i),psi_scan);
    AF2(:,i) = array.compute_array_factor(v,v2(:,i),psi_scan);
%     AF(:,i) = array.compute_array_factor(v,v1(:,i)+v2(:,i),psi_scan);
%    
%    AF(:,i) = array.compute_array_factor(v,x(:,i),psi_scan);
end

fig6=figure;
%plot(u_scan,(abs(AF)./max(abs(AF),[],1)));
hold on;
plot(u_scan,20*log10(abs(AF1)./max(abs(AF1),[],1)),'k-');
plot(u_scan,20*log10(abs(AF2)./max(abs(AF2),[],1)),'k-');
plot(u1*[1 1],[-20 0],'k--');
plot(u2*[1 1],[-20 0],'k--');
xlabel('u');

annotation(fig6,'doublearrow',[0.477430555555556 0.560763888888889],...
    [0.252086419753086 0.251543209876543]);

text(0,-18,'$\delta_u$','FontSize',11);
grid off

ylabel('Array Factor [linear]');
ylim([-20 0]);
%legend(arrayfun(@(x) sprintf('d = %.1f \\delta_u',x),spacing,'UniformOutput',false));
utils.setPlotStyle(gca,{'widescreen'});
utils.exportPlot(fig6,[prefix '6']);

%% Figure 7 - Array Tapers

% Initialize tapers
N = 11;
w_u = utils.makeTaper(N,'uniform');
w_c = utils.makeTaper(N,'cosine');
w_h = utils.makeTaper(N,'hann');
w_hh = utils.makeTaper(N,'hamming');

w = [w_u;w_c;w_h;w_hh];

fig7a=figure;
hold on;
for idx_plot = 1:4
    set(gca,'LineStyleOrderIndex',idx_plot);
    this_width = 2 - .25*idx_plot;
    plot(0:N-1,w(idx_plot,:),'linewidth',this_width);
end
legend('uniform','cosine','hann','hamming');
xlabel('i');
ylabel('$a_i$');
grid off;
ylim([0 1.1]);
xlim([-0.5 10.5])
utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig7a,[prefix '7a']);

% Generate beampatterns

% Over-sampling factor; for smoothness
osf = 100;
u = linspace(-1,1,osf*N);

% Take the fourier transform, and normalize peak response
W = fftshift(fft(w,osf*N,2),2);
W = W ./ max(abs(W),[],2);

fig7b=figure;hold on;
for idx_plot = 1:4
    set(gca,'LineStyleOrderIndex',idx_plot);
    this_width = 2 - .25*idx_plot;
    plot(u,20*log10(abs(W(idx_plot,:))),'linewidth',this_width);
end
legend('uniform','cosine','hann','hamming');
xlim([0 1]);
ylim([-80 0]);
xlabel('$u = \cos(\theta)$')
ylabel('$|G(\theta)|^2$ [dB]');
grid off;

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(fig7b,[prefix '7b']);

%% Figure 9/11 - Beamscan Images
% Set up array
N = 25;
d_lam = .5;
d_psi = .89/((N-1)*d_lam);
du = sin(d_psi);
v = array.make_steering_vector(d_lam,N);

% Set up source signals
spacing = 1;
u1 = du*spacing/2;
u2 = -u1;

psi1 = asin(u1);
psi2 = asin(u2);

v1 = v(psi1); % N x numel(u1)
v2 = v(psi2); % N x numel(u2)

% Generate snapshots
M = 30;
s1 = sqrt(.5)*(randn(1,M)+1i*randn(1,M));
s2 = sqrt(.5)*(randn(1,M)+1i*randn(1,M));
x1 = v1.*s1;
x2 = v2.*s2;

% Add noise
snr_db = 10;
snr = 10.^(snr_db/10);
n = sqrt(1./(snr*2))*(randn(size(x1))+1i*randn(size(x1)));
x = x1 + x2 + n;


% Generate Beamscan images
[P,psi_vec] = array.beamscan(x,v,pi/2,1001);
P_mvdr = array.beamscan_mvdr(x,v,pi/2,1001);

% Scale outputs
P = P/max(P(:));
P_mvdr = P_mvdr/max(P_mvdr(:));

fig9=figure;hold on;
plot(sin(psi_vec),10*log10(P),'LineWidth',1.5,'DisplayName','Beamscan');
set(gca,'LineStyleOrderIndex',2);
plot(sin(psi_vec),10*log10(P_mvdr),'LineWidth',1.25,'DisplayName','MVDR');
%h1=plot(u1*[1 1],[-60 0],'k--');
%h2=plot(u2*[1 1],[-60 0],'k--');
%utils.excludeFromLegend([h1 h2]);
ylim([-30 0]);
xlabel('u');
ylabel('P [dB]');
legend('Location','NorthEast');
utils.setPlotStyle(gca,{'widescreen'});
utils.exportPlot(fig9,[prefix '9']);

P_music = array.music(x,v,2,pi/2,1001);
P_music = P_music./max(abs(P_music(:)));
set(gca,'LineStyleOrderIndex',3);
plot(sin(psi_vec),10*log10(abs(P_music)),'DisplayName','MUSIC');
utils.setPlotStyle(gca,{'widescreen'});
utils.exportPlot(fig9,[prefix '11']);

%% Figure 10 - Beamscan Example Images

fig10=ex8_1;
utils.exportPlot(fig10,[prefix '10']);

%% Figure 12 - CRLB Plot
if force_recalc
    
xi_dB = -20:.5:0;
xi_lin = 10.^(xi_dB/10);
psi = 5*pi/180;
M = 100;
N = 11;
[v,v_dot] = array.make_steering_vector(.5,N);

% Compute CRLB
C_psi = zeros(numel(xi_dB),1);
C_psi_stoch = zeros(numel(xi_dB),1);

for idx_xi = 1:numel(xi_dB)
    C_psi(idx_xi) = array.crlb_det(xi_lin(idx_xi),1,psi,M,v,v_dot);
    C_psi_stoch(idx_xi) = array.crlb_stochastic(xi_lin(idx_xi),1,psi,M,v,v_dot);
end
C_th = (180/pi)^2*C_psi;
RMSE_th = sqrt((C_th));
C_th_stoch =(180/pi)^2*C_psi_stoch;
RMSE_th_stoch = sqrt((C_th_stoch));

figcrlb=figure;
semilogy(xi_dB,RMSE_th,'k--','DisplayName','Det. CRLB');
hold on;
plot(xi_dB,RMSE_th_stoch,'k-.','DisplayName','Stochastic CRLB');
legend('Location','NorthEast');

nMC = 1000;
s = exp(1i*2*pi*rand(1,M,nMC));%sqrt(1/2)*(randn(1,M,nMC)+1i*randn(1,M,nMC));
x0 = v(psi).*s;
n = sqrt(1/2)*(randn(N,M,nMC)+1i*randn(N,M,nMC));
RMSE_th_beam = zeros(numel(xi_dB),1);
RMSE_th_mvdr = zeros(numel(xi_dB),1);
RMSE_th_music = zeros(numel(xi_dB),1);
fprintf('Executing array DF monte carlo trial...\n\t');
tStart = tic;
for idx_xi = 1:numel(xi_dB)
    fprintf('.');
    sn2 = 1./xi_lin(idx_xi);
    x = x0 + sqrt(sn2)*n;
    
    this_err_beamscan = zeros(1,nMC);
    this_err_mvdr = zeros(1,nMC);
    this_err_music = zeros(1,nMC);

    for idxMC=1:nMC
        % Compute beamscan image
        [P,psi_vec] = array.beamscan(x(:,:,idxMC),v,pi/2,2001);
        [~,idx_pk] = max(abs(P));
        this_err_beamscan(idxMC) = abs(psi_vec(idx_pk)-psi);
    
        % Compute beamscan image
        [P_mvdr,psi_vec] = array.beamscan_mvdr(x(:,:,idxMC),v,pi/2,2001);
        [~,idx_pk_mvdr] = max(abs(P_mvdr));
        this_err_mvdr(idxMC) = abs(psi_vec(idx_pk_mvdr)-psi);
    
        % Compute beamscan image
        [P_music,psi_vec] = array.music(x(:,:,idxMC),v,1,pi/2,2001);
        [~,idx_pk_music] = max(abs(P_music));
        this_err_music(idxMC) = abs(psi_vec(idx_pk_music)-psi);
    end
        
    % Average Results
    RMSE_th_beam(idx_xi) = (180/pi)*sqrt(sum(this_err_beamscan.^2)/nMC);
    RMSE_th_mvdr(idx_xi) = (180/pi)*sqrt(sum(this_err_mvdr.^2)/nMC);
    RMSE_th_music(idx_xi) = (180/pi)*sqrt(sum(this_err_music.^2)/nMC);
        
end

fprintf('done.\n');
tElapsed = toc(tStart);
hrs = floor(tElapsed/3600);
mins = floor((tElapsed-hrs*3600)/60);
secs = tElapsed - 60*mins - 3600*hrs;
fprintf('Elapsed time: %d hours, %d minutes, %.2f seconds\n',hrs,mins,secs);

figure(figcrlb);
set(gca,'ColorOrderIndex',2);
plot(xi_dB,RMSE_th_beam,'DisplayName','Beamscan','LineWidth',2);
plot(xi_dB,RMSE_th_mvdr,'DisplayName','MVDR','LineWidth',1.5);
plot(xi_dB,RMSE_th_music,'DisplayName','MUSIC','LineWidth',1);
xlabel('SNR [dB]');
ylabel('RMSE [deg]');
grid on;

utils.setPlotStyle(gca,{'widescreen','tight'});
utils.exportPlot(figcrlb,[prefix '12']);
end

%% Figure 13 - Example 8.2
if force_recalc

fig13=ex8_2;
utils.exportPlot(gcf,[prefix '13']);

end

%% Cleanup

% Restore plot settings
utils.resetPlotSettings;