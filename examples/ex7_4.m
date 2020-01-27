function fig = ex7_4(fig_in)
% fig=ex7_4(fig_in)
%
% Executes Example 7.4.  If the input fig_in is provided, then the curve is
% added to the existing plot, otherwise a new figure is generated.
%
% INPUTS
%   fig_in      figure handle to pre-existing figure, on which new data
%               will be plotted
%
% OUTPUTS
%   fig         figure handle to generated graphic
%
% Nicholas O'Donoughue
% 1 July 2019

if nargin<1 || isempty(fig_in)
    fig = figure;
end
% Parameter Definition
Pt = 35e-3; % W
Gt = 34; % dBi
B_s = 31.3e6; % Hz -- unused
f0 = 35e9; % Hz
lambda = utils.constants.c/f0; % m
Lt = 0; % dB
ht = 500; % m -- unused
Lr = 2; % dB
F_n = 4; % dB
B_n = 40e6; % Hz
T_s = 10e-6; % s
hr = 500; % m

% Angle of Arrival
th_true = 10;
psi_true = th_true*pi/180;

% Compute the SNR, as a function of range
R = [1e3:1e3:10e3, 20e3:10e3:100e3];
S = (10*log10(Pt) + Gt - Lt) + 20*log10(lambda./(4*pi*R)) - Lr;
N = 10*log10(utils.constants.boltzmann*utils.constants.T0*B_n)+F_n;
snr_db = S-N;

% Compute the number of time samples available
M = 1 + floor(T_s*B_n);

% Narrow Baseline Interferometer
d_lam = .5;
rmse_max = 1;
crlb_max = (rmse_max*pi/180)^2;
snr_min_eff = (1./(2*M*crlb_max)).*(1./(2*pi*d_lam*cos(psi_true))).^2;
snr_min_lin = 2*snr_min_eff;
R2 = 1e3:10:10e3;
S2 = (10*log10(Pt) + Gt - Lt) + 20*log10(lambda./(4*pi*R2)) - Lr;
N2 = 10*log10(utils.constants.boltzmann*utils.constants.T0*B_n)+F_n;
snr_db2 = S2-N2;
[~,idx] = min(abs(snr_db2 - 10*log10(snr_min_lin)));
R_interf_narrowBaseline = R2(idx);

% Wide Baseline Interferometer
d_lam = 2;
rmse_max = 1;
crlb_max = (rmse_max*pi/180)^2;
snr_min_eff = (1./(2*M*crlb_max)).*(1./(2*pi*d_lam*cos(psi_true))).^2;
snr_min_lin = 2*snr_min_eff;
R2 = [1e3:100:40e3];
S2 = (10*log10(Pt) + Gt - Lt) + 20*log10(lambda./(4*pi*R2)) - Lr;
N2 = 10*log10(utils.constants.boltzmann*utils.constants.T0*B_n)+F_n;
snr_db2 = S2-N2;
[~,idx] = min(abs(snr_db2 - 10*log10(snr_min_lin)));
R_interf_wideBaseline = R2(idx);

% Compute RMSE
d_lam = .5;
crlb_psi_interf = aoa.interf_crlb(snr_db,snr_db,M,d_lam,psi_true);
rmse_th_interf = (180/pi)*sqrt(crlb_psi_interf);

d_lam = 2;
crlb_psi_interf = aoa.interf_crlb(snr_db,snr_db,M,d_lam,psi_true);
rmse_th_interf2 = (180/pi)*sqrt(crlb_psi_interf);

% Add to plot
hold on;
loglog(R/1e3,rmse_th_interf,'-o','DisplayName','Interferometer (d=\lambda/2)');
loglog(R/1e3,rmse_th_interf2,'-x','DisplayName','Interferometer (d=2\lambda)');
legend('location','NorthWest');
