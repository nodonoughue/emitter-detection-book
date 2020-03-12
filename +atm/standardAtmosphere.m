function atmStruct = standardAtmosphere(alt)
% atmStruct = standardAtmosphere(alt)
%
% Implement the standard atmosphere recommended
% by ITU-R P.835-6 (12/2017).
%
% Requires an input altitude.  If an array of altitudes
% are defined, the function is run once per input and the
% array of outputs are returned.
%
% Inputs:
%   alt: Altitude [m], can be scalar or N-dimensional matrix
%
% Outputs:
%   T:   Temperature [K]
%   P:   Total pressure [hPa]
%   rho: Water Vapor Density [g/m^3]
%   e:   Wator Vapor Partial Pressure [hPa]
%
% Currently using the mean annual standard temperature,
% as opposed to the seasonal and latitude-specific models
% also provided in ITU-R P.835-6 (12/2017).
%
% Nicholas O'Donoughue
% 1 July 2019

if nargin < 1 || isempty(alt)
    alt = 0; % Just above sea level
end

%% Use recursion to handle array inputs
if numel(alt) > 1
    % Use array fun to call standardAtmosphere on each altitude
    structs = arrayfun(@(x) atm.standardAtmosphere(x),alt);
    
    % Convert response from an array of structs to a struct with array
    % fields
    dims = size(alt);
    atmStruct = struct('alt',alt,...
                       'P',reshape([structs(:).P],dims),...
                       'T',reshape([structs(:).T],dims),...
                       'e',reshape([structs(:).e],dims),...
                       'R',reshape([structs(:).R],dims),...
                       'M',reshape([structs(:).M],dims));
    return;
end

% At this point, we can assume alt is a scalar

%% Compute Geopotential Height from Geometric Height
% Equation 1a
hp = 6356.766*(alt/1e3) ./ (6356.766+(alt/1e3)); % [km']
h = 6356.766*hp ./ (6356.766-hp);

%% Temperature (K)
% Equation (2a)-(2g)
if hp <= 11
    T = 288.15-6.5*hp;
elseif hp<= 20
    T = 216.65;
elseif hp <= 32
    T = 216.65 + (hp-20);
elseif hp <= 47
    T = 228.65 + 2.8 * (hp-32);
elseif hp <= 51
    T = 270.65;
elseif hp <= 71
    T = 270.65 - 2.8*(hp-51);
elseif hp <= 84.852
    T = 214.65-2*(hp-71);
elseif h <= 91 % Region 2
    T = 186.8673;
elseif h <= 100
    T = 263.1905-76.3232*(1-((h-91)/19.9429)^2)^.5;
else
    error('Standard temperature not defined above 100 km altitude.');
end
    
%% Pressure (hPa)
if hp <= 11
    P = 1013.25*(288.15/(288.15-6.5*hp))^(-34.1623/6.5);
elseif hp <= 20
    P = 226.3226*exp(-34.1634*(hp-11)/216.65);
elseif hp <= 32
    P = 54.74980*(216.65/(216.65+(hp-20)))^34.1632;
elseif hp <= 47
    P = 8.680422*(228.65/(228.65+2.8*(hp-32)))^(34.1632/2.8);
elseif hp <= 51
    P = 1.109106*exp(-34.1632*(hp-47)/170.65);
elseif hp <= 71
    P = 0.6694167*(270.65/(270.65-2.8*(hp-51)))^(-34.1632/2.8);
elseif hp <= 84.852
    P = 0.03956649*(214.65/(214.65-2*(hp-71)))^(-34.1632/2);
elseif h <= 100
    a0 = 95.571899;
    a1 = -4.011801;
    a2 = 6.424731e-2;
    a3 = -4.789660e-4;
    a4 = 1.340543e-6;
    P = exp(a0+a1*h+a2*h^2+a3*h^3+a4*h^4);
else
    error('Standard pressure not defined above 100 km altitude.');
end

%% Water Vapor Pressure/Density
h0 = 2;     % Scale height [km]
rho0 = 7.5; % Standard ground-level density [g/m^3]
rho = rho0*exp(-h/h0); % [g/m^3]

e = rho*T/216.67; % [hPa]

% Adjust water vaper pressure for constant mixing above a certain altitude
 mixingRatio = max(2e-6,e./P); % See text following equation 8
 e = P.*mixingRatio;
 
 %% Wrap it in a struct
 atmStruct = struct('alt',h,...
                    'P',P,...
                    'T',T,...
                    'e',e,...
                    'R',0,...
                    'M',0);