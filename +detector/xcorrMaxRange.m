function R = xcorrMaxRange(PFA,PD,Tcorr,Tp,Bn,Bs,f0,ht,hr,SNR0,includeAtmLoss,atmStruct)
% R = xcorrMaxRange(PFA,PD,Tcorr,Tp,Bn,Bs,f0,ht,hr,SNR0,includeAtmLoss,
%                                                       atmStruct)
%
% Compute the maximum range for a square law detector, as specified by the
% PD, PFA, and number of samples (M).  The link is described by the carrier
% frequency (f0), and transmit/receive antenna heights (ht and hr), and the
% transmitter and receiver are specified by the SNR in the absence of path
% loss (SNR0).  If specified, the atmospheric struct is passed onto the
% path loss model.
%
% Inputs:
%
%   PFA             Probability of False Alarm [0-1]
%   PD              Probability of Detection [0-1]
%   Tcorr           Correlation time [sec]
%   Tp              Pulse Duration [sec]
%   Bn              Noise bandwidth [Hz]
%   Bs              Signal Bandwidth [Hz]
%   f0              Carrier Frequency [Hz]
%   ht              Transmitter height [m]
%   hr              Receiver height [m]
%   SNR0            SNR in the absence of path loss [dB]
%   includeAtmLoss  Boolean flag indicating whether atmospheric loss should
%                   be considered in path loss [Default = false]
%   atmStruct       (Optional) struct containing atmospheric parameters;
%                   can be called from atm.standardAtmosphere().
%
% Outputs:
%   R               Maximum range at which the specified PD/PFA condition
%                   can be met [m]
%
% Nicholas O'Donoughue
% 1 July 2019

if nargin < 6 || isempty(includeAtmLoss)
    includeAtmLoss = false;
    atmStruct = [];
end

% Find the required SNR Threshold
SNRmin = detector.xcorrMinSNR(PFA,PD,Tcorr,Tp,Bn,Bs);

R = zeros(size(SNRmin));

for ii=1:numel(SNRmin)
    thisSNRmin = SNRmin(ii);
    if numel(SNR0) > 1
        thisSNR0 = SNR0(ii);
    else
        thisSNR0 = SNR0;
    end
    
    
    % Find the acceptable propagation loss
    Lprop_max = thisSNR0-thisSNRmin;
    
    % Set up error function
    Lprop_fun = @(R) prop.pathLoss(R,f0,ht,hr,includeAtmLoss,atmStruct);
    err_fun = @(R) Lprop_fun(R)-Lprop_max;
    
    % Set up initial search point
    thisR = 1e3;
    err = err_fun(thisR);
    
    % Optimization Parameters
    err_tol = .01;      % SNR error tolerance [dB]
    maxIter = 1000;     % Maximum number of iterations
    iter = 0;           % Iteration counter
    
    % Perform the optimization
    while iter < maxIter && abs(err)>err_tol
        % Compute derivative
        dR = 1; % 1 meter
        y1 = err_fun(thisR+dR);
        y0 = err_fun(thisR-dR);
        df = (y1-y0)/(2*dR);
        
        % Newton Method step
        thisR = thisR - err/df;
        err = err_fun(thisR);
        
        % Iteration count
        iter = iter+1;
    end
    R(ii) = thisR;
end