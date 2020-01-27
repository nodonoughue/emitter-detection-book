% Example approach to analyze an interferometer
%
% 25 July 2019
% Nicholas O'Donoughue

% Generate the Signals
th_true = 45;
d_lam = .5;
psi_true = th_true*pi/180;
phi = 2*pi*d_lam*sin(psi_true);
alpha = 1;

% Set up the parameter sweep
M_vec = [10,100];       % Number of temporal samples at each antenna 
                        % test point
snr_db_vec = -10:.2:20; % signal to noise ratio
num_mc = 1e4;           % number of monte carlo trials at each parameter 
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
    
    % Generate Signals
    s1 = sqrt(.5)*(randn(M,this_num_mc)+1i*randn(M,this_num_mc));
    s2 = alpha*s1*exp(1i*phi);
    
    % Generate Noise
    noise_base1 = sqrt(1/2)*(randn(M,this_num_mc)+1i*randn(M,this_num_mc));
    noise_base2 = sqrt(1/2)*(randn(M,this_num_mc)+1i*randn(M,this_num_mc));

    % Loop over SNR levels
    for idx_snr = 1:numel(snr_db_vec)
        if mod(idx_snr,10)==0
            fprintf('.');
        end
        
        % Compute noise power, scale base noise
        snr_db = snr_db_vec(idx_snr);
        noise_pwr = 10.^(-snr_db/10);
        
        % Generate noisy signals
        x1 = s1 + sqrt(noise_pwr)*noise_base1;
        x2 = s2 + sqrt(noise_pwr)*noise_base2;
        
        % Compute the estimate for each Monte Carlo trial
        psi_est = zeros(1,this_num_mc);
        for idx_mc = 1:this_num_mc
            psi_est(idx_mc) = aoa.interf_df(x1(:,idx_mc),x2(:,idx_mc),d_lam);
        end
        
        % Compute RMS Error
        rmse_psi(idx_M,idx_snr) = sqrt(sum(abs((psi_est-psi_true)).^2)/this_num_mc);
        
        % Compute CRLB for RMS Error
        crlb_psi(idx_M,idx_snr) = aoa.interf_crlb(snr_db,snr_db+20*log10(alpha),M,d_lam,psi_true);
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
title('Interferometer DF Performance');
legend('Location','SouthWest');

text(5,7,'M=10','FontSize',10);
text(5,2,'M=100','FontSize',10);