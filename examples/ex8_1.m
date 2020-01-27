function fig=ex8_1()
% fig=ex8_1()
%
% Executes Example 8.1, relying on the sample data in examples/ex8_1.mat,
% and generates one figure.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig     figure handle to generated graphic
%
% Nicholas O'Donoughue
% 1 July 2019

% Load sample data
load('examples/ex8_1.mat','x','D','N','d_lam');
%   x       noisy data vector (M x N)
%   D       number of sources
%   N       number of array elements
%   M       number of time snapshots
%   d_lam   array spacing

% Construct array steering vector
v = array.make_steering_vector(d_lam,N);

% Call Beamformer
[P,psi_vec] = array.beamscan(x,v,pi/2,1001);

[peak_val,peak_idx] = findpeaks(P,psi_vec);
[~,sort_idx] = sort(peak_val,'descend');
peak_vals = peak_val(sort_idx(1:D));
psi_soln = peak_idx(sort_idx(1:D));
th_soln = 180*psi_soln/pi;

% Call MVDR Beamformer
[P_mvdr,psi_vec] = array.beamscan_mvdr(x,v,pi/2,1001);
[peak_val_mvdr,peak_idx] = findpeaks(P_mvdr,psi_vec);
[~,sort_idx] = sort(peak_val_mvdr,'descend');
peak_vals_mvdr = peak_val_mvdr(sort_idx(1:D));
psi_soln_mvdr = peak_idx(sort_idx(1:D));
th_soln_mvdr = 180*psi_soln_mvdr/pi;

% Plot
th_vec = 180*psi_vec/pi;
fig=figure;
plot(th_vec,10*log10(abs(P)),'LineWidth',1.5,'DisplayName','Beamscan');hold on;
set(gca,'ColorOrderIndex',3);
plot(th_vec,10*log10(abs(P_mvdr)),'DisplayName','MVDR');
plot(th_soln,10*log10(peak_vals),'kv','DisplayName','Beamscan Soln.','MarkerSize',6);
plot(th_soln_mvdr,10*log10(peak_vals_mvdr),'k^','DisplayName','MVDR Soln.','MarkerSize',6);
%utils.excludeFromLegend(h0(2:end));
% utils.excludeFromLegend(h1(2:end));
xlabel('$\theta$ [deg]');
ylabel('P [dB]');
ylim([-145 -100]);
legend('Location','NorthWest');

utils.setPlotStyle(gca,{'widescreen','tight'});
