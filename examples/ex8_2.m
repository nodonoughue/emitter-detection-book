function fig=ex8_2()
% fig = ex8_2()
%
% Executes Example 8.2 and generates one figure
%
% INPUTS
%   none
%
% OUTPUTS
%   fig     figure handle to generated graphic
%
% Nicholas O'Donoughue
% 1 July 2019

% Define transmitter values
Pt = 100;
Gt = 0;
Bs = 500e3;
f0 = 950e6;
Lt = 3;
ht = 5e3;
theta = 5;%[0 30 60]+rand(1,3)/10;
psi = theta*pi/180;

% Define receiver values
d_lam = .5;
N = 5;
Gr = 0;
Lr = 2;
Fn = 4;
Bn = 100e6;
hr = 10;
T = 1e-6;

% Compute the number of snapshots
t_samp = 1/(2*f0);
M = floor(T/t_samp);

% Generate the array factor
[v,v_dot] = array.make_steering_vector(d_lam,N);

% Generate received signal power and SNR at multiple ranges
R = 100e3:20e3:500e3;
Lprop = prop.pathLoss(R,f0,ht,hr,false);
Pr_dB = 10*log10(Pt) + Gt + Gr - Lt - Lr - Lprop;
N_dB = 10*log10(utils.constants.kT*Bn) + Fn;
xi_dB = Pr_dB-N_dB;
xi_lin = 10.^(xi_dB/10);
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

% Compute MC Experiment
nMC = 1000;
s = sqrt(1/2)*(randn(M,nMC)+1i*randn(M,nMC));
n = sqrt(1/2)*(randn(N,M,nMC)+1i*randn(N,M,nMC));
RMSE_th_beam = zeros(numel(R),numel(psi));
RMSE_th_mvdr = zeros(numel(R),numel(psi));
RMSE_th_music = zeros(numel(R),numel(psi));
fprintf('Executing array DF monte carlo trial...\n\t');
tStart = tic;
for idxR = 1:numel(R)
    fprintf('.');
    thisS = s*sqrt(10.^(xi_dB(idxR)/10));
    
    for idxPsi = 1:numel(psi)
        
        thisX = v(psi(idxPsi)).*reshape(thisS,1,M,nMC) + n;
        
        this_err_beamscan = zeros(1,nMC);
        this_err_mvdr = zeros(1,nMC);
        this_err_music = zeros(1,nMC);
        for idxMC=1:nMC
            % Compute beamscan image
            [P,psi_vec] = array.beamscan(thisX(:,:,idxMC),v,pi/2,2001);
            [~,idx_pk] = max(abs(P));
            this_err_beamscan(idxMC) = abs(psi_vec(idx_pk)-psi(idxPsi));
            
            % Compute beamscan image
            [P_mvdr,psi_vec] = array.beamscan_mvdr(thisX(:,:,idxMC),v,pi/2,2001);
            [~,idx_pk] = max(abs(P_mvdr));
            this_err_mvdr(idxMC) = abs(psi_vec(idx_pk)-psi(idxPsi));
            
            % Compute beamscan image
            [P_music,psi_vec] = array.music(thisX(:,:,idxMC),v,1,pi/2,2001);
            [~,idx_pk] = max(abs(P_music));
            this_err_music(idxMC) = abs(psi_vec(idx_pk)-psi(idxPsi));
%             if idxMC==1
%                 figure;
%                 plot(psi_vec,10*log10(abs(P)),psi_vec,10*log10(abs(P_mvdr)),psi_vec,10*log10(abs(P_music)));
%                 legend('Beamscan','MVDR','MUSIC');
%                 xlabel('$\psi$');
%                 keyboard;
%             end
        end
        
        % Average Results
        RMSE_th_beam(idxR,idxPsi) = (180/pi)*sqrt(sum(this_err_beamscan.^2)/nMC);
        RMSE_th_mvdr(idxR,idxPsi) = (180/pi)*sqrt(sum(this_err_mvdr.^2)/nMC);
        RMSE_th_music(idxR,idxPsi) = (180/pi)*sqrt(sum(this_err_music.^2)/nMC);
        
    end
end  
fprintf('done.\n');
tElapsed = toc(tStart);
hrs = floor(tElapsed/3600);
mins = floor((tElapsed-hrs*3600)/60);
secs = tElapsed - 60*mins - 3600*hrs;
fprintf('Elapsed time: %d hours, %d minutes, %.2f seconds\n',hrs,mins,secs);

% Plot results
fig=figure;
loglog(R/1e3,RMSE_th,'k--','DisplayName','CRLB (det.)');
hold on;
loglog(R/1e3,RMSE_th_stoch,'k-.','DisplayName','CRLB (stoch.)');
set(gca,'ColorOrderIndex',2);
plot(R/1e3,RMSE_th_beam,'LineWidth',2,'DisplayName','Beamscan');
plot(R/1e3,RMSE_th_mvdr,'LineWidth',1.5,'DisplayName','MVDR');
plot(R/1e3,RMSE_th_music,'LineWidth',1,'DisplayName','MUSIC');
xlabel('Range [km]');
ylabel('RMSE [deg]');
grid on;
legend('Location','NorthWest')
utils.setPlotStyle(gca,{'widescreen','tight'});