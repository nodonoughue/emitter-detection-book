function make_itu_rec_figures
% Generates figures to compare with ITU-R P.676-12
%
% Can be used to ensure atmospheric loss tables and calculations are
% reasonably accurate.
%
% Nicholas O'Donoughue
% 1 July 2019

%% Figure 1 - Specific attenuation
f_ghz = 0:1:1000;
[f_ctr_o,f_ctr_w] = atm.getSpectralLines(f_ghz(1)*1e9,f_ghz(end)*1e9);
f = sort([f_ghz * 1e9,f_ctr_o,f_ctr_w]);
f_ghz = f/1e9; % recompute f_ghz to include spectral lines
atmStruct = atm.standardAtmosphere(0);

[gamma_ox,gamma_h2o] = atm.gasLossCoeff(f,atmStruct.P,atmStruct.e,atmStruct.T);

figure;
semilogy(f_ghz,gamma_ox,'b-','DisplayName','Dry');
hold on;
plot(f_ghz,gamma_ox+gamma_h2o,'r-','DisplayName','Standard');
xlabel('Frequency (GHz)');
ylabel('Specific Attenuation (dB/km)');
grid on;
legend('Location','NorthWest');
title('Replication of ITU-R P.676-12, Figure 1');
xlim([0 1e3]);

%% Figure 2
f_ghz = 50:.1:70;
[f_ctr_o,f_ctr_w] = atm.getSpectralLines(f_ghz(1)*1e9,f_ghz(end)*1e9);
f = sort([f_ghz * 1e9,f_ctr_o,f_ctr_w]);
f_ghz = f/1e9; % recompute f_ghz to include spectral lines

alts = (0:5:20)*1e3;
atmStruct = atm.standardAtmosphere(alts);

[gamma_ox,gamma_h2o] = atm.gasLossCoeff(f(:),atmStruct.P,atmStruct.e,atmStruct.T);
gamma = gamma_ox + gamma_h2o;

figure;
semilogy(f_ghz,gamma,'-');
legend(arrayfun(@(x) sprintf('%d km',x/1e3),alts,'UniformOutput',false));
xlabel('Frequency (GHz)');
ylabel('Specific Attenuation (dB/km)');
grid on;
title('Replication of ITU-R P.676-12, Figure 2');

%% Figure 4
f_ghz = 0:1:1e3;
[f_ctr_o,f_ctr_w] = atm.getSpectralLines(f_ghz(1)*1e9,f_ghz(end)*1e9);
f = sort([f_ghz * 1e9,f_ctr_o,f_ctr_w]);
f_ghz = f/1e9; % recompute f_ghz to include spectral lines

[l,lo,lw] = atm.calcZenithLoss(f,0,0);

figure;
semilogy(f_ghz,lo,'b-','DisplayName','Dry');
hold on;
semilogy(f_ghz,l,'r-','DisplayName','Standard');
xlim([0 1e3]);
xlabel('Frequency (Ghz)');
ylabel('Zenith Attenuation (dB)');
grid on;
legend('Location','NorthWest');
title('Replication of ITU-R P.676-12, Figure 4');

%% Figure 10

f_ghz = 1:350;
[f_ctr_o,f_ctr_w] = atm.getSpectralLines(f_ghz(1)*1e9,f_ghz(end)*1e9);
f = sort([f_ghz * 1e9,f_ctr_o,f_ctr_w]);
f_ghz = f/1e9; % recompute f_ghz to include spectral lines

atmStruct = atm.standardAtmosphere(0);

[gamma_ox,gamma_h2o] = atm.gasLossCoeff(f,atmStruct.P,atmStruct.e,atmStruct.T);

figure;
loglog(f_ghz,gamma_ox,'b-','DisplayName','Dry');
hold on;
plot(f_ghz,gamma_h2o,'r-','DisplayName','Water Vapour');
plot(f_ghz,gamma_ox+gamma_h2o,'k-','DisplayName','Total');
xlabel('Frequency (GHz)');
ylabel('Specific Attenuation (dB/km)');
grid on;
legend('Location','NorthWest');
title('Replication of ITU-R P.676-12, Figure 10');
xlim([0 350]);

%% Figure 11

f_ghz = 1:350;
[f_ctr_o,f_ctr_w] = atm.getSpectralLines(f_ghz(1)*1e9,f_ghz(end)*1e9);
f = sort([f_ghz * 1e9,f_ctr_o,f_ctr_w]);
f_ghz = f/1e9; % recompute f_ghz to include spectral lines

[l,lo,lw] = atm.calcZenithLoss(f,0,0);

figure;
loglog(f_ghz,lo,'b-','DisplayName','Dry');
hold on;
loglog(f_ghz,lw,'r-','DisplayName','Water vapour');
semilogy(f_ghz,l,'k-','DisplayName','Total');
xlim([0 350]);
xlabel('Frequency (Ghz)');
ylabel('Zenith Attenuation (dB)');
grid on;
legend('Location','NorthWest');
title('Replication of ITU-R P.676-12, Figure 11');

%% Figure 12

f_ghz = 50:.01:70;
[f_ctr_o,f_ctr_w] = atm.getSpectralLines(f_ghz(1)*1e9,f_ghz(end)*1e9);
f = sort([f_ghz * 1e9,f_ctr_o,f_ctr_w]);
f_ghz = f/1e9; % recompute f_ghz to include spectral lines

alts = (0:5:20)*1e3;

l = atm.calcZenithLoss(f(:),alts,0);

figure;
loglog(f_ghz,l,'-');%'DisplayName',sprintf('%d km',alts(idx_alt)/1e3));
xlabel('Frequency (GHz)');
ylabel('Zenith Attenuation (dB)');
grid on;
legend(arrayfun(@(x) sprintf('%d km',x/1e3),alts,'UniformOutput',false));
title('Replication of ITU-R P.676-12, Figure 12');
