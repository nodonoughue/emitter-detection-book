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

% Coeffs for kh
a = [-5.3398,-0.35351,-0.23789,-0.94158];
b = [-0.10008,1.26970,0.86036,0.64552];
c = [1.13098,0.454,0.15354,0.16817];
m = -0.18961;
ck = 0.71147;

log_kh = sum(a(:).*exp(-((log10(f/1e9)-b(:))./c(:)).^2))+m*log10(f/1e9)+ck;
kh = 10.^(log_kh);


% Coeffs for kv
a = [-3.80595,-3.44965,-0.39902,0.50167];
b = [0.56934,-0.22911,0.73042,1.07319];
c = [0.81061,0.51059,0.11899,0.27195];
m = -0.16398;
ck = 0.63297;

log_kv = sum(a(:).*exp(-((log10(f/1e9)-b(:))./c(:)).^2))+m*log10(f/1e9)+ck;
kv = 10.^(log_kv);


% Coeffs for ah
a = [-0.14318,0.29591,0.32177,-5.37610,16.1721];
b = [1.82442,0.77564,0.63773,-0.96230,-3.29980];
c = [-0.55187,0.19822,0.13164,1.47828,3.43990];
m = 0.67849;
ca=-1.95537;

ah = sum(a(:).*exp(-((log10(f/1e9)-b(:))./c(:)).^2))+m*log10(f/1e9)+ca;

% Coeffs for av
a = [-0.07771,0.56727,-0.20238,-48.2991,48.5833];
b = [2.33840,0.95545,1.14520,0.791669,0.791459];
c = [-0.76284,0.54039,0.26809,0.116226,0.116479];
m = -0.053739;
ca = 0.83433;

av = sum(a(:).*exp(-((log10(f/1e9)-b(:))./c(:)).^2))+m*log10(f/1e9)+ca;

%% Account for Polarization and Elevation Angles
k = .5*(kh + kv + (kh-kv)*cos(el_angle).^2*cos(2*pol_tilt));
a = (kh.*ah+kv.*av+(kh.*ah-kv.*av).*cos(el_angle).^2.*cos(2*pol_tilt))./(2*k);

gamma = k.*rainfall_rate.^a;

