function fig=ex5_2
% fig=ex5_2()
%
% Executes Example 5.2.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle to generated graphic
%
% Nicholas O'Donoughue
% 1 July 2019

Thop = [1e-3,1e-4,1e-5]';
Bs = 1e5:1e5:1e7;
Bhop = 4e9;

Tdwell = 1./Bs;

Nscan = Thop./Tdwell;
Br = max(Bs,Bhop./Nscan);

fig=figure;
colors = get(0,'DefaultAxesColorOrder');
set(0,'DefaultAxesColorOrder',colors(1,:));
loglog(Bs/1e6,Br/1e6);
xlabel('Frequency Resolution ($\delta_f$) [MHz]');
ylabel('Receiver Bandwidth ($B_r$) [MHz]');
legend(arrayfun(@(x) sprintf('T_{hop}=%.2f ms',x*1e3),Thop,'UniformOutput',false));

utils.setPlotStyle(gca,{'widescreen','tight'});