function L = pathLoss(R,f0,ht,hr,includeAtmLoss,atmStruct)
% L = pathLoss(R,f0,ht,hr,includeAtmLoss,atmStruct)
%
% Computes the propagation loss according to a piece-wise
% model where free space is used at close range, and two-ray
% is used at long range.  The cross-over range between the two is
% the Fresnel Zone.
%
% Inputs:
%   R                   Range of link (m)
%   f0                  Carrier frequency (Hz)
%   ht                  Transmitter height (m)
%   hr                  Receiver height (m)
%   includeAtmLoss      Boolean flag.  If true (default) then
%                       atmospheric absorption is modeled and
%                       returned for a standard atmosphere.
%   atmStruct           Atmospheric loss parameter struct, must match the
%                       format expected by calcAtmLoss.
%
% Outputs:
%   L                   Path Loss [dB]
%
% Nicholas O'Donoughue
% 1 July 2019

if nargin < 5 || isempty(includeAtmLoss)
    includeAtmLoss = true;
end

% Find the fresnel zone distance
fz = prop.fresnelZone(f0,ht,hr);

% Compute free space path loss - w/out atmospherics
fspl = prop.freeSpacePathLoss(R,f0,false);
tworay = prop.twoRayPathLoss(R,f0,ht,hr,false);

% Combine the free space and two ray path loss calculations,
% using binary singleton expansion to handle non-uniform
% parameter sizes, so long as all non-singleton dimension match,
% this will succeed.
fsplMask = bsxfun(@le,R,fz);
L = bsxfun(@times,fspl,fsplMask) + ...
    bsxfun(@times,tworay,~fsplMask);
       
if includeAtmLoss
    if nargin < 6 || isempty(atmStruct)
        atmStruct = atm.standardAtmosphere(.5*(ht+hr));
    end
    atmLoss = atm.calcAtmLoss(f0,R,0,0,atmStruct);
    L = bsxfun(@plus,L,atmLoss);
end

