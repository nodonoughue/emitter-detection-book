function [rng_bistatic, rr_bistatic] = measurement(x_tx, x_rx, x_tgt, v_tx, v_rx, v_tgt)
% 
%
% Computed bistatic range and range rate (Doppler) measurements for a 
% set of bistatic transmitter/receiver pairs.
%
% INPUTS:
%   x_tx             nDim x nTx array of transmitter coordinates
%   x_rx             nDim x nRx array of receiver coordinates
%   x_tgt            nDim x nTgt array of target coordinates
%   v_tx             nDim x nTx array of transmitter velocities
%   v_rx             nDim x nRx array of receiver velocities
%   v_tgt            nDim x nTgt array of target velocities
%
% OUTPUTS:
%   rng_bistatic     nTx x nRx x nTgt matrix of bistatic range
%   rr_bistatic      nTx x nRx x nTgt matrix of bistatic range-rate
%
% Nicholas O'Donoughue
% 13 March 2020

do_doppler = nargin >= 4;
    
% Parse inputs
[nDim,nTx] = size(x_tx);
[~,nRx] = size(x_rx);
[~,nTgt] = size(x_tgt);

%% Bistatic Range
% Compute distance from each source position to each sensor
dtx = reshape(x_tgt,nDim,1,1,nTgt) - reshape(x_tx,nDim,nTx);
drx = reshape(x_tgt,nDim,1,1,nTgt) - reshape(x_rx,nDim,1,nRx);

Rtx = sqrt(sum(abs(dtx).^2,1)); % 1 x nTx x 1 x nTgt
Rrx = sqrt(sum(abs(drx).^2,1)); % 1 x 1 x nRx x nTgt

% Assemble Bistatic Range
rng_bistatic = reshape(Rtx + Rrx, nTx, nRx, nTgt);

%% Bistatic Range Rate
if do_doppler
    if isempty(v_tx)
        % Default is no motion
        v_tx = zeros(size(x_tx));
    end

    if isempty(v_rx)
        % Default is no motion
        v_rx = zeros(size(x_rx));
    end

    if isempty(v_tgt)
        % Default is no motion
        v_tgt = zeros(size(x_tgt));
    end

    % Compute range rate from range and velocity
    vtx = reshape(v_tgt,nDim,1,1,nTgt) - reshape(v_tx,nDim,nTx);
    vrx = reshape(v_tgt,nDim,1,1,nTgt) - reshape(v_rx,nDim,1,nRx);

    rrtx = reshape(-sum(vtx.*dtx./Rtx,1),nTx,1,nTgt);
    rrrx = reshape(-sum(vrx.*drx./Rrx,1),1,nRx,nTgt);

    % Assemble bistatic range rate
    rr_bistatic = rrtx + rrrx;
else
    rr_bistatic = zeros(size(rng_bistatic));
end