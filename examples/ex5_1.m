function fig = ex5_1()
% fig=ex5_1()
%
% Executes Example 5.1.
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

% Scan Variables
Thop = [10e-3,1e-3,1e-4]; % Target signal hopping period
Bhop = 200e6; % Target signal hopping bandwidth
Bs = 5e6;     % Target signal transmit bandwidth
Br = 5e6;    % Target signal receive bandwidth

Tdet = Thop*Br/Bhop;

fs = Br;
M = floor(Tdet*fs);

% Detection Curve
xi_db_vec = -15:.1:10;
xi_lin_vec = 10.^(xi_db_vec/10);
PFA = 1e-6;

[XI,MM] = ndgrid(xi_lin_vec,M);
if use_stat_toolbox
    eta = chi2inv(1-PFA,2*MM);
    PD = 1-ncx2cdf(eta,2*MM,2*MM.*XI);
else
    eta = utils.chi2inv(1-PFA,2*MM);
    PD = 1-ncx2cdf(eta,2*MM,2*MM.*XI);
end

colors = get(0,'DefaultAxesColorOrder');
set(0,'DefaultAxesColorOrder',colors(1,:));
fig=figure;
plot(xi_db_vec,PD);
xlabel('SNR [dB]');ylabel('$P_D$');
legend(arrayfun(@(x) sprintf('T_{hop}=%.1f ms',x),Thop*1e3,'UniformOutput',false));
set(0,'DefaultAxesColorOrder',colors);
utils.setPlotStyle(gca,{'widescreen','tight'});
