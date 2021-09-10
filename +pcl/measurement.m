function rho = measurement(x_tx, x_rx, x_tgt, v_tx, v_rx, v_tgt, ref_idx)
% 
%
% Computed bistatic range and range rate (Doppler) measurements for a 
% set of bistatic transmitter/receiver pairs.
%
% If the velocity terms are all blank, then range only measurements are
% generated.  Otherwise range-rate measurements are generated, and returned
% stacked vertically after the range measurements.
%
% INPUTS:
%   x_tx             nDim x nTx array of transmitter coordinates
%   x_rx             nDim x nRx array of receiver coordinates
%   x_tgt            nDim x nTgt array of target coordinates
%   v_tx             nDim x nTx array of transmitter velocities
%   v_rx             nDim x nRx array of receiver velocities
%   v_tgt            nDim x nTgt array of target velocities
%   ref_idx          Matrix of tx/rx pairing indices (in the order
%                    they're used in C).  If ignored, then all pairwise
%                    measurements are used (nTx x nRx)
%
% OUTPUTS:
%   rho              nMsmt x nTgt matrix of bistatic range (and range-rate)
%                    measurements.
%
% Nicholas O'Donoughue
% 13 March 2020

do_doppler = nargin >= 4 && (~isempty(v_tx) || ~isempty(v_rx) || ~isempty(v_tgt))

if nargin < 7 || ~exist('ref_idx','var')
    ref_idx = [];
end

if do_doppler && isempty(v_tx)
    v_tx = zeros(size(x_tx));
end

if do_doppler && isempty(v_rx)
    v_rx = zeros(size(x_rx));
end

if do_doppler && isempty(v_tgt)
    v_tgt = zeros(size(x_tgt));
end

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
else
    rr_bistatic = [];
end


%% Parse Tx/Rx pairing indices
if isempty(ref_idx) || strcmpi(ref_idx,'full')==0
    % Nothing (or full set of pairs) specified, do all pairs
    rng_bistatic = reshape(Rtx + Rrx, nTx*nRx, nTgt);
    
    if do_doppler
        rr_bistatic = reshape(rrtx + rrrx, nTx*nRx, nTgt);
    end
elseif isscalar(ref_idx)
    % Scalar specified, just keep that tx-rx pair
    rng_bistatic = reshape(Rtx(1,ref_idx,1,:) + ...
                           Rrx(1,1,ref_idx,:), 1, nTgt);
    
    if do_doppler
        rr_bistatic = reshape(rrtx(1,ref_idx,1,:) + ...
                              rrrx(1,1,ref_idx,:), 1, nTgt);
    end
else
    % Matrix; it must have two rows.  Each column is a tx/rx pair setting
    % transmitters indices on the first row, receiver on the second
    assert(size(ref_idx,1)==2,'Error parsing reference index; incorrect size.')
    assert(ismatrix(ref_idx),'Error parsing reference index; incorrect shape.');
    
    tx_idx = ref_idx(1,:);
    rx_idx = ref_idx(2,:);
    
    rng_bistatic = reshape(squeeze(Rtx(:,tx_idx,:,:) + Rrx(:,:,rx_idx,:)),...
                           numel(tx_idx), nTgt);
    
    if do_doppler
        rr_bistatic = reshape(squeeze(rrtx(:,tx_idx,:,:) + rrrx(:,:,rx_idx,:)),...
                              numel(tx_idx), nTgt);
    end
    
end

rho = cat(1,rng_bistatic, rr_bistatic);