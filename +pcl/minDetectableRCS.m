function min_rcs = minDetectableRCS(x_tx, x_rx, x_tgt, tx_fom, rx_fom, snr_min, ref_idx)
% function min_rcs = minDetectableRCS(x_tx, x_rx, x_tgt, tx_fom, rx_fom, snr_min, ref_idx)
%
% Computes the min detectable signal, given one or more bistatic geometries
% using PCL.
%
%
% Inputs:
%   x_tx                Transmitter positions [m] - nDim x nTx
%   x_rx                Receiver positions [m] - nDim x nRx
%   x_tgt               Target positions to test [m] - nDim x nTgt
%   tx_fom              Transmit FOM [W*Hz*m^3] (EIRP*BW*lam^3/L) - 1 x nTx
%   rx_fom              Receiver FOM [S] (Gr*T_int/(P_noise*Loss)
%   snr_min             Required SNR [dB] for detection
%   ref_idx             Matrix of tx/rx pairing indices (in the order
%                       they're used in C).  If ignored, then all pairwise
%                       measurements are used (nTx x nRx)
% Outputs:
%   min_rcs             Min Detect RCS [dBsm]
%
% 13 September 2021
% Nicholas O'Donoughue


[nDim,nTx] = size(x_tx);
[nDim2,nRx] = size(x_rx);
[nDim3,nTgt] = size(x_tgt);
assert(nDim==nDim2 && nDim==nDim3,'All input positions must have the same number of physical dimensions.');
assert(isscalar(tx_fom) || numel(tx_fom)==nTx,'Number of transmit FOMs must be scalar, or equal to the number of transmit positions.');
assert(isscalar(rx_fom) || numel(rx_fom)==nRx,'Number of receive FOMs must be scalar, or equal to the number of receive positions.');

%% Bistatic Range
% Compute distance from each source position to each sensor
dtx = reshape(x_tgt,nDim,1,1,nTgt) - reshape(x_tx,nDim,nTx);
drx = reshape(x_tgt,nDim,1,1,nTgt) - reshape(x_rx,nDim,1,nRx);

Rtx = sqrt(sum(abs(dtx).^2,1)); % 1 x nTx x 1 x nTgt
Rrx = sqrt(sum(abs(drx).^2,1)); % 1 x 1 x nRx x nTgt

%% Figure of Merits
fom = reshape(tx_fom,1,nTx,1,1) .* reshape(rx_fom,1,1,nRx,1);

%% Parse Tx/Rx pairing indices
if isempty(ref_idx) || strcmpi(ref_idx,'full')==0
    % Nothing (or full set of pairs) specified, do all pairs
    Rtx = reshape(repmat(Rtx,1,1,nRx,1),[nDim,nTx*nRx,nTgt]);
    Rrx = reshape(repmat(Rrx,1,nTx,1,1),[nDim,nTx*nRx,nTgt]);
    fom = reshape(fom,[1,nTx*nRx,1]);
    
elseif isscalar(ref_idx)
    % Scalar specified, just keep that tx-rx pair
    Rtx = reshape(Rtx(:,ref_idx,1,:),[nDim,1,nTgt]);
    Rrx = reshape(Rrx(:,1,ref_idx,:),[nDim,1,nTgt]);
    fom = squeeze(fom(1,ref_idx,ref_idx,1));
    
else
    % Matrix; it must have two rows.  Each column is a tx/rx pair setting
    % transmitters indices on the first row, receiver on the second
    assert(size(ref_idx,1)==2,'Error parsing reference index; incorrect size.')
    assert(ismatrix(ref_idx),'Error parsing reference index; incorrect shape.');
    
    tx_idx = ref_idx(1,:);
    rx_idx = ref_idx(2,:);
    
    Rtx = reshape(Rtx(:,tx_idx,1,:),[nDim,numel(tx_idx),nTgt]);
    Rrx = reshape(Rrx(:,1,rx_idx,:),[nDim,numel(rx_idx),nTgt]);
    fom = reshape(fom(1,tx_idx,rx_idx,1),[1,numel(rx_idx),1]);
    
end

%% Min Detectable RCS
min_rcs_lin = 10.^(snr_min/10) * (4*pi)^3 .* Rtx.^2 .* Rrx.^2 ./ fom;
min_rcs = 10*log10(min_rcs_lin);
