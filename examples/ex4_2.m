function fig=ex4_2()
% fig=ex4_2()
%
% Executes Example 4.2.
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
Pt = 10*log10(500); % tower / user
f0 = 16e9;
ht = 6096;
Lt = 3;
Bs = 200e6;
Tp = 20e-6;
Gt = [30,10,0];

% Receive Side
hr = 20;
Gr = 0;
Lr = 3;
NF = 4;
Bn = 500e6;
Tcorr = 1e-4;


% Compute xi0
N0 = utils.constants.boltzmann*utils.constants.T0*10^(NF/10);
N = 10*log10(N0*Bn);
xi0 = Pt+Gt + Gr - Lt - Lr - N;

% Adjust xi0 to account for partial pulse
xi0 = xi0 + 10*log10(min(1,Tp/Tcorr));

% Compute Prop Loss
Rvec = 10e3:10e3:1000e3;
Lprop = prop.pathLoss(Rvec,f0,ht,hr,false);

xi = xi0(:) - Lprop;

% Compute thresholds
PD=.8;
PFA=1e-6;
M = fix(Tcorr*Bn);
xi_ed = detector.squareLawMinSNR(PFA,PD,M);
xi_xc = detector.xcorrMinSNR(PFA,PD,Tcorr,Tp,Bn,Bs);

% Compute Max Range
R_ml_ed = detector.squareLawMaxRange(PFA,PD,M,f0,ht,hr,xi0(1),false,[]);
R_nearSL_ed = detector.squareLawMaxRange(PFA,PD,M,f0,ht,hr,xi0(2),false,[]);
R_farSL_ed = detector.squareLawMaxRange(PFA,PD,M,f0,ht,hr,xi0(3),false,[]);
R_ml_xc = detector.xcorrMaxRange(PFA,PD,Tcorr,Tp,Bn,Bs,f0,ht,hr,xi0(1),false,[]);
R_nearSL_xc = detector.xcorrMaxRange(PFA,PD,Tcorr,Tp,Bn,Bs,f0,ht,hr,xi0(2),false,[]);
R_farSL_xc = detector.xcorrMaxRange(PFA,PD,Tcorr,Tp,Bn,Bs,f0,ht,hr,xi0(3),false,[]);

fprintf('Mainlobe detection\n');
fprintf('\tusing Energy Detector: %.2f km\n',R_ml_ed/1e3);
fprintf('\tusing Cross Correlator: %.2f km\n',R_ml_xc/1e3);
fprintf('Near sidelobe detection\n');
fprintf('\tusing Energy Detector: %.2f km\n',R_nearSL_ed/1e3);
fprintf('\tusing Cross Correlator: %.2f km\n',R_nearSL_xc/1e3);
fprintf('Far sidelobe detection\n');
fprintf('\tusing Energy Detector: %.2f km\n',R_farSL_ed/1e3);
fprintf('\tusing Cross Correlator: %.2f km\n',R_farSL_xc/1e3);

% Plot results
fig=figure;
colors = get(0,'DefaultAxesColorOrder');
set(gca,'ColorOrder',colors([1,3,4],:));
hold on;
legendCell = {'\xi (ML)','\xi (Near SL)','\xi (Far SL)'};
for i=1:numel(Gt)
    set(gca,'ColorOrderIndex',i);
    plot(Rvec/1e3,xi(i,:),'DisplayName',legendCell{i});
end
plot(Rvec([1 end])/1e3,xi_ed*[1 1],'k--','DisplayName','\xi_{min} (ED)');
plot(Rvec([1 end])/1e3,xi_xc*[1 1],'k-.','DisplayName','\xi_{min} (XC)');
legend;
grid on;
xlabel('Range [km]');
ylabel('$\xi$ [dB]');
set(gca,'xscale','log');
utils.setPlotStyle(gca,{'widescreen','tight'});
