function fig=ex4_1()
% fig=ex4_1()
%
% Executes Example 4.1.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle to generated graphic
%
% Nicholas O'Donoughue
% 1 July 2019

% Transmit Side
ERP = [10 -7]; % tower / user
f0 = [850e6 1900e6];
ht = [60 2];
Lt = 0;
Bs = 5e6;
Tp = 20e-3;

% Receive Side
hr = 1e3;
Gr = 0;
Lr = 3;
NF = 5;
Bn = 100e6;
Tcorr = 1e-4;


% Compute xi0
N0 = utils.constants.boltzmann*utils.constants.T0*10^(NF/10);
N = 10*log10(N0*Bn);
xi0 = ERP + Gr - Lt - Lr - N;

% Compute Prop Loss
Rvec = 1e3:1e3:1000e3;
Lprop_twr = prop.pathLoss(Rvec,f0(1),ht(1),hr,false);
Lprop_user= prop.pathLoss(Rvec,f0(2),ht(2),hr,false);

xi_twr = xi0(1) - Lprop_twr;
xi_user = xi0(2) -Lprop_user;

% Compute thresholds
PD=.8;
PFA=1e-6;
M = fix(Tcorr*Bn);
xi_ed = detector.squareLawMinSNR(PFA,PD,M);
xi_xc = detector.xcorrMinSNR(PFA,PD,Tcorr,Tp,Bn,Bs);

% Compute Max Range
R_twr_ed = detector.squareLawMaxRange(PFA,PD,M,f0(1),ht(1),hr,xi0(1),false,[]);
R_user_ed = detector.squareLawMaxRange(PFA,PD,M,f0(2),ht(2),hr,xi0(2),false,[]);
R_twr_xc = detector.xcorrMaxRange(PFA,PD,Tcorr,Tp,Bn,Bs,f0(1),ht(1),hr,xi0(1),false,[]);
R_user_xc = detector.xcorrMaxRange(PFA,PD,Tcorr,Tp,Bn,Bs,f0(2),ht(2),hr,xi0(2),false,[]);

fprintf('Detection range of tower\n');
fprintf('\tusing Energy Detector: %.2f km\n',R_twr_ed/1e3);
fprintf('\tusing Cross Correlator: %.2f km\n',R_twr_xc/1e3);
fprintf('Detection range of handset user\n');
fprintf('\tusing Energy Detector: %.2f km\n',R_user_ed/1e3);
fprintf('\tusing Cross Correlator: %.2f km\n',R_user_xc/1e3);

% Plot results
fig=figure;
semilogx(Rvec/1e3,xi_twr,'DisplayName','\xi (Tower)');hold on;
set(gca,'ColorOrderIndex',4);
plot(Rvec/1e3,xi_user,'DisplayName','\xi (Handset)');
plot(Rvec([1 end])/1e3,xi_ed*[1 1],'k--','DisplayName','\xi_{min} (ED)');
plot(Rvec([1 end])/1e3,xi_xc*[1 1],'k-.','DisplayName','\xi_{min} (XC)');
legend;
grid on;
xlabel('Range [km]');
ylabel('$\xi$ [dB]');

utils.setPlotStyle(gca,{'widescreen','tight'});
