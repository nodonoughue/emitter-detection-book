function figs = book2_ex3_1()
% fig=book2_ex3_1()
%
% Executes Example 3.1 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   figs         array of figure handles
%
% Nicholas O'Donoughue
% 28 May 2021

%% Initialize Measurement Covariance Matrix
C = diag([1, 3, 2, 3, 5].^2); % Square to convert from std dev to variance
c_max = 5^2 + 3^2;

%% Generate Common Reference Sets
C_first = utils.resampleCovMtx(C, 1);
C_last = utils.resampleCovMtx(C, 5);

figa = figure();
imagesc(C_first);
caxis([0 c_max]);
colorbar;
title('Ref Index = 1');

figb = figure();
imagesc(C_last);
caxis([0 c_max]);
colorbar;
title('Ref Index = 5');

%% Generate Full Measurement Set
C_full = utils.resampleCovMtx(C, 'full');

figc = figure();
imagesc(C_full);
caxis([0 c_max]);
colorbar;
title('Full Measurement Set');

%% Output Figures
figs = [figa, figb, figc];