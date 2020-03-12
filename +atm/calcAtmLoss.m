function L = calcAtmLoss(f,Dg,Dr,Dc,atmStruct,pol_angle,el_angle)
% L = calcAtmLoss(f,Dg,Dr,Dc,atmStruct,pol_angle,el_angle)
%
% INPUTS:
%   f               - Frequency [Hz]
%   Dg              - Path length for gas loss [m]
%   Dr              - Path length for rain loss [m]
%   Dc              - Path length for cloud loss [m]
%   atmStruct       - Struct with atmospheric parameters
%                       'P' - Dry air pressure [hPa]
%                       'e' - Water vapor partial pressure [hPa]
%                       'T' - Temperature [K]
%                       'R' - Rainfall rate [mm/hr]
%                       'M' - Cloud density [g/m^3]
%   pol_angle       - Polarization angle [radians], 0 for Horizontal,
%                     pi/2 for Vertical, between 0 and pi for slant.
%                     [default = 0]
%   el_angle        - Elevation angle of the path under test [default = 0]
%   
%
% OUTPUT:
%   L               - Loss along path due to atmospheric absorption
%
%
% Ref:
%   ITU-R P.676-11(09/2016) Attenuation by atmospheric gases
%   ITU-R P.840-6 (09/2013) Attenuation due to clouds and fog
%   ITU-R P.838-3 (03/2005) Specific attenuation model for rain for use in
%   prediction methods
%
%
% Nicholas O'Donoughue
% 1 July 2019

%% Check optional inputs
if nargin < 7 || isempty(el_angle)
    el_angle = 0;
end

if nargin < 6 || isempty(pol_angle)
    pol_angle = 0;
end

if nargin < 5 || isempty(atmStruct)
    % Default atmosphere is the standard atmosphere at sea level, with no
    % fog/clouds or rain.
    atmStruct = atm.standardAtmosphere(0);
end

if nargin < 4 || isempty(Dc)
    Dc = 0;
end

if nargin < 3 || isempty(Dr)
    Dr = 0;
end

if isempty(Dg)
    Dg = 0;
end

%% Compute loss coefficients
if Dg > 0
    [gammaOx,gammaW] = atm.gasLossCoeff(f,atmStruct.P,atmStruct.e,atmStruct.T);
    gamma_g = gammaOx+gammaW;
else
    gamma_g = 0;
end

if Dr > 0 && atmStruct.R > 0
    gamma_r = atm.rainLossCoeff(f,pol_angle,el_angle,atmStruct.R);
else
    gamma_r = 0;
end

if Dc > 0 && atmStruct.M > 0
    gamma_c = atm.fogLossCoeff(f,atmStruct.M,atmStruct.T);
else
    gamma_c = 0;
end

%% Compute loss componenets
L_g = gamma_g .* Dg/1e3;
L_r = gamma_r .* Dr/1e3;
L_c = gamma_c .* Dc/1e3;

L = L_g + L_r + L_c;