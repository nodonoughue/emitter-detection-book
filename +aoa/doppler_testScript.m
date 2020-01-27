% Example analysis of a Doppler DF receiver
%
% 25 July 2019
% Nicholas O'Donoughue

% Generate the Signals
th_true = 45;
psi_true = th_true*pi/180;
A = 1;
phi0 = 2*pi*rand();
f = 1e9;
ts = 1/(5*f);

% Doppler antenna parameters
psi_0 = 0;          % initial pointing angle [radians]
c=3e8;              % speed of light [m/s]
lam = c/f;          % wavelength [m]
R = lam/2;          % Doppler radius (chosen to be half-wavelength) [m]
psi_res = .0001;    % Desired DF resolution [rad]

% Set up the parameter sweep
M_vec = [10,100,1000];  % Number of temporal samples at each antenna 
                        % test point
snr_db_vec = -10:2:20;  % signal to noise ratio
num_mc = 1e6;           % number of monte carlo trials at each parameter 
                        % setting

% Set up output variables
rmse_psi = zeros(numel(M_vec),numel(snr_db_vec));
crlb_psi = zeros(numel(M_vec),numel(snr_db_vec));

% Loop over parameters
fprintf('Executing Doppler Monte Carlo sweep...\r\n');
for idx_M = 1:numel(M_vec)
    M = M_vec(idx_M);
    this_num_mc = num_mc / M;
    fprintf('\t M=%d',M);
    
    % Reference signal
    t_vec = ts*(0:M-1);
    r0 = A*exp(1i*phi0)*exp(1i*2*pi*f*t_vec);
    
    % Doppler signal
    fr = 1/(ts*M); % Ensure a single cycle during M
    x0 = A*exp(1i*phi0)*exp(1i*2*pi*f*t_vec).*exp(1i*2*pi*f*R/c*cos(2*pi*fr*t_vec-psi_true));
    
    % Generate noise signal
    noise_base_r = sqrt(sum(abs(r0).^2)/(M*2))*(randn(M,this_num_mc)+1i*randn(1,this_num_mc));
    noise_base_x = sqrt(sum(abs(r0).^2)/(M*2))*(randn(M,this_num_mc)+1i*randn(1,this_num_mc));

    % Loop over SNR levels
    for idx_snr = 1:numel(snr_db_vec)
         if mod(idx_snr,10)==0
            fprintf('.');
         end
        
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
        
        % Compute RMS Error
        rmse_psi(idx_M,idx_snr) = sqrt(sum(abs((psi_est-psi_true)).^2)/this_num_mc);
        
        % Compute CRLB for RMS Error
        crlb_psi(idx_M,idx_snr) = aoa.doppler_crlb(snr_db,M,A,ts,f,R,fr,psi_true);
    end
fprintf('done.\n');
end


figure;
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
