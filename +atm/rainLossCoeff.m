function gamma = rainLossCoeff(f,pol_tilt,el_angle,rainfall_rate)
% gamma = rainLossCoeff(f,pol_tilt,el_angle,rainfall_rate)
%
% Computes the rain loss coefficient given a frequency, polarization,
% elevation angle, and rainfall rate, according to ITU-R P.838-3, 2005.
%
% Inputs:
%   f               - Propagation Frequency [Hz]
%   pol_tilt        - Polarization angle [radians], 0 = Horizontal and 
%                     pi/2 is Vertical.  Slanted polarizations will have
%                     a value 0 and pi.
%   el_angle        - Propagation path elevation angle [radians]
%   rainfall_rate   - Rainfall rate [mm/hr]
%
% Outputs:
%   gamma           - Loss coefficient [dB/km] caused by rain.
%
%
% Nicholas O'Donoughue
% 1 July 2019

% Find first unused dimension
n_dims_used = max([numel(size(f)),numel(size(pol_tilt)),numel(size(el_angle)),numel(size(rainfall_rate))]);
coeff_dim = n_dims_used+1;

% Coeffs for kh
a = reshape([-5.3398,-0.35351,-0.23789,-0.94158],[ones(1,n_dims_used),4]);
b = reshape([-0.10008,1.26970,0.86036,0.64552],[ones(1,n_dims_used),4]);
c = reshape([1.13098,0.454,0.15354,0.16817],[ones(1,n_dims_used),4]);
m = -0.18961;
ck = 0.71147;

log_kh = sum(a.*exp(-((log10(f/1e9)-b)./c).^2),coeff_dim)+m*log10(f/1e9)+ck;
kh = 10.^(log_kh);


% Coeffs for kv
a = reshape([-3.80595,-3.44965,-0.39902,0.50167],[ones(1,n_dims_used),4]);
b = reshape([0.56934,-0.22911,0.73042,1.07319],[ones(1,n_dims_used),4]);
c = reshape([0.81061,0.51059,0.11899,0.27195],[ones(1,n_dims_used),4]);
m = -0.16398;
ck = 0.63297;

log_kv = sum(a.*exp(-((log10(f/1e9)-b)./c).^2),coeff_dim)+m*log10(f/1e9)+ck;
kv = 10.^(log_kv);


% Coeffs for ah
a = reshape([-0.14318,0.29591,0.32177,-5.37610,16.1721],[ones(1,n_dims_used),5]);
b = reshape([1.82442,0.77564,0.63773,-0.96230,-3.29980],[ones(1,n_dims_used),5]);
c = reshape([-0.55187,0.19822,0.13164,1.47828,3.43990],[ones(1,n_dims_used),5]);
m = 0.67849;
ca=-1.95537;

ah = sum(a.*exp(-((log10(f/1e9)-b)./c).^2),coeff_dim)+m*log10(f/1e9)+ca;

% Coeffs for av
a = reshape([-0.07771,0.56727,-0.20238,-48.2991,48.5833],[ones(1,n_dims_used),5]);
b = reshape([2.33840,0.95545,1.14520,0.791669,0.791459],[ones(1,n_dims_used),5]);
c = reshape([-0.76284,0.54039,0.26809,0.116226,0.116479],[ones(1,n_dims_used),5]);
m = -0.053739;
ca = 0.83433;

av = sum(a.*exp(-((log10(f/1e9)-b)./c).^2),coeff_dim)+m*log10(f/1e9)+ca;

%% Account for Polarization and Elevation Angles
k = .5*(kh + kv + (kh-kv).*cos(el_angle).^2.*cos(2*pol_tilt));
a = (kh.*ah+kv.*av+(kh.*ah-kv.*av).*cos(el_angle).^2.*cos(2*pol_tilt))./(2*k);

gamma = k.*rainfall_rate.^a;

