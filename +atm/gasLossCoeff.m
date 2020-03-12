function [gammaOx,gammaW] = gasLossCoeff(f,p,e,T)
% [gammaOx,gammaW] = gasLossCoeff(f,P,e,T)
%
% Inputs:
%   f       - Propagation Frequencies [Hz]
%   P       - Dry Air Pressure [hPa]
%   e       - Water Vapor Partial Pressure [hPa]
%   T       - Temperature [K]
%
% Outputs:
%   gammaOx - Gas loss coefficient due to oxygen [dB/km]
%   gammaW  - Gas loss coefficient due to water vapor [dB/km]
%
% Implement the atmospheric loss coefficients from Annex 1 of ITU-R
% P.676-11 (12/2017)
%
% If array inputs are specified, then array results are given for
% alphaO and alphaW.
%
% 1 July 2019
% Nicholas O'Donoughue

if all(f(:)==0 | isnan(f(:)))
    gammaOx=0;
    gammaW=0;
    return;
end

% Determine largest dimension in use
n_dims_in_use = max([numel(size(f)),numel(size(p)),numel(size(e)),numel(size(T))]);
spec_dim = n_dims_in_use + 1;

% Read in the spectroscopic tables (Tables 1 and 2 of Annex 1)
% All table data will be Nx1 for the N spectroscopic lines of
% each table.
[fox,a1,a2,a3,a4,a5,a6] = atm.makeSpectroscopicTableOxygen();
[fw,b1,b2,b3,b4,b5,b6] = atm.makeSpectroscopicTableWater();

% Shift spectroscopic data to use the first non-singleton dimension across
% all inputs
reshape_dims = [ones(1,n_dims_in_use),numel(fox)];
fox = reshape(fox,reshape_dims);
a1 = reshape(a1,reshape_dims);
a2 = reshape(a2,reshape_dims);
a3 = reshape(a3,reshape_dims);
a4 = reshape(a4,reshape_dims);
a5 = reshape(a5,reshape_dims);
a6 = reshape(a6,reshape_dims);

reshape_dims = [ones(1,n_dims_in_use),numel(fw)];
fw = reshape(fw,reshape_dims);
b1 = reshape(b1,reshape_dims);
b2 = reshape(b2,reshape_dims);
b3 = reshape(b3,reshape_dims);
b4 = reshape(b4,reshape_dims);
b5 = reshape(b5,reshape_dims);
b6 = reshape(b6,reshape_dims);

%% Compute the dry continuum due to pressure-induced Nitrogen absorption and the Debye spectrum (eq 8)
f0 = f/1e9; % Convert freq from Hz to GHz
th = 300./T;
d = 5.6e-4*(p+e).*th.^.8;
ND = f0.*p.*th.^2.*(6.14e-5./(d.*(1+(f0./d).^2))+(1.4e-12*p.*th.^1.5)./(1+1.9e-5*f0.^1.5));

%% Compute the strength of the i-th water/o2 vapor line (eq 3)
Sox = (a1.*1e-7.*p.*th.^3).*exp(a2.*(1-th));
Sw = (b1.*1e-1.*e.*th.^3.5).*exp(b2.*(1-th));

%% Compute the line shape factor for each
%  Correction factor due to interference effects in oxygen lines (eq 7)
dox = (a5+a6.*th).*1e-4.*(p+e).*th.^.8;
dw = 0;

% spectroscopic line width (eq 6a)
dfox = (a3*1e-4).*(p.*th.^(.8-a4)+1.1*e.*th);
dfw = (b3*1e-4).*(p.*th.^(b4)+b5.*e.*th.^b6);

% modify spectroscopic line width to account for Zeeman splitting of oxygen
% lines and Doppler broadening of water vapour lines (eq 6b)
dfox_sq = dfox.^2+2.25e-6;
dfox = sqrt(dfox_sq);
dfw = .535.*dfw+sqrt(.217*dfw.^2+(2.1316e-12.*fw.^2)./th);

% Compute line shape factor
delta_fox = fox-f0;
sum_fox = fox+f0;
delta_fw = fw-f0;
sum_fw = fw+f0;
Fox = f0./fox.*( ((dfox-dox.*delta_fox)./(delta_fox.^2+dfox.^2)) + ((dfox-dox.*sum_fox)./(sum_fox.^2+dfox_sq)));
Fw = f0./fw.*( ((dfw-dw.*delta_fw)./(delta_fw.^2+dfw.^2)) + ((dfw-dw.*sum_fw)./(sum_fw.^2+dfw.^2)));

%% Compute complex refractivities
Nox = sum(Sox.*Fox,spec_dim)+ND;
Nw = sum(Sw.*Fw,spec_dim);

gammaOx = .1820*f0.*Nox;
gammaW = .1820*f0.*Nw;

%% Handle all freqs < 1 GHz
lowFreq = f0 < 1;
gammaOx(lowFreq) = 0;
gammaW(lowFreq) = 0;
