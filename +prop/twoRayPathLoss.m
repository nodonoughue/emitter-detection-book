function L = twoRayPathLoss(R,f0,ht,hr,includeAtmLoss)
% L = twoRayPathLoss(R,f0,ht,hr,includeAtmLoss)
%
% Computes the two-ray path loss according to
%     L = 10*log10(R^4/(h_t^2*h_r^2))
%
% This model is generally deemed appropriate for low altitude transmitters
% and receivers, with a flat Earth model.
%
% Inputs:
%
%   R                   Range (meters)
%   f0                  Carrier frequency (Hz)
%   ht                  Height of the transmitter (m)
%   hr                  Height of the receiver (m)
%   includeAtmLoss      Boolean flag indicating whether
%                       atmospheric loss should be included
%                       (Default = True)
%
% Outputs:
%   L                   Loss [dB]
%
% If the input includeAtmLoss is True (default),
% then a call is also made to the atmLoss function, and
% the total loss is returned.
%
% All non-scalar inputs must have the same dimensions
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 5 || isempty(includeAtmLoss)
    includeAtmLoss = true;
end

tworay = 10*log10(R.^4./(ht.^2.*hr.^2));

if includeAtmLoss
    atmStruct = atm.standardAtmosphere(.5*(ht+hr));
    atmLoss = atm.calcAtmLoss(f0,R,0,0,atmStruct);
else
    atmLoss = 0;
end

L = tworay+atmLoss;