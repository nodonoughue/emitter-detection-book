function L = freeSpacePathLoss(R,f0,includeAtmLoss,atmStruct)
% L = freeSpacePathLoss(R,f0,includeAtmLoss)
%
% Computes the free space path loss according to:
%     L = 20*log10(4*pi*R/lambda)
%
% Inputs:
%   R                   Range [m]
%   f0                  Carrier frequency [Hz]
%   includeAtmLoss      Boolean flag indicating whether
%                       atmospheric loss should be included
%                       (Default = True)
%   atmStruct           Atmospheric loss parameter struct, must match the
%                       format expected by calcAtmLoss.
%
% Outputs:
%   L                   Path loss [dB]
%
% If the input includeAtmLoss is True (default),
% then a call is also made to the atmLoss function, and
% the total loss is returned.
%
% All non-scalar inputs must have the same dimensions
%
% Nicholas O'Donoughue
% 1 July 2019

if nargin < 3 || isempty(includeAtmLoss)
    includeAtmLoss = true;
end

% Convert from frequency to wavelength [m]
lam = utils.constants.c./f0;

% Equation (B.1)
fspl = 20*log10(4*pi*R./lam);

% Add atmospheric loss, if called for
if includeAtmLoss
    if nargin < 4 || isempty(atmStruct)
        atmStruct = atm.standardAtmosphere(.5*(ht+hr));
    end
    atmLoss = atm.calcAtmLoss(f0,R,0,0,atmStruct);
else
    atmLoss = 0;
end

L = fspl+atmLoss;