function gamma = fogLossCoeff(f,M,T)
% gamma = fogLossCoeff(f,M,T)
%
% Implement the absorption loss coefficient due to clouds and fog, as a
% function of the frequency, cloud density, and temperature, according to
% ITU-R P.840-7 (2017).
%
% Inputs:
%   f       - Propagation Frequencies [Hz]
%   M       - Cloud/fog density [g/m^3]
%   T       - Atmospheric temperature [K]
%
% Outputs:
%   gamma   - Loss coefficient [dB/km]
%
% Nicholas O'Donoughue
% 1 July 2019

if nargin < 3 || isempty(T)
    atmStruct = atm.standardAtmosphere;
    T = atmStruct.T;
end

%% Cloud Liquid Water Specific Attenuation Coefficient
theta = 300./T;
e0=77.66+103.3*(theta-1);
e1=0.0671*e0;
e2=3.52;

fp = 20.20-146*(theta-1)+316*(theta-1).^2;
fs = 39.8*fp;

e_prime= (e0-e1)./(1+((f/1e9)./fp).^2)+(e1-e2)./(1+((f/1e9)./fs).^2)+e2;
e_prime_prime= (f/1e9).*(e0-e1)./(fp.*(1+(f/1e9./fp).^2))+((f/1e9).*(e1-e2)./(fs.*(1+((f/1e9)./fs).^2)));

eta = (2+e_prime)./(e_prime_prime);
Kl = .819*(f/1e9)./(e_prime_prime.*(1+eta.^2));

%% Cloud attenuation
gamma = Kl.*M;