function J = jacobian(x_tx, x_rx, x_tgt, v_tx, v_rx, v_tgt, ref_idx)
% J = jacobian(x_tx, x_rx, x_source, v_tx, v_rx, v_source, ref_idx)
%
% Returns the Jacobian matrix for PCL bistatic range and range-rate
% measurements of a source at x_source (nDim x nSource).
%
% TODO: Incorporate range-rate measurements
% TODO: Incorporate AOA processing
%
% INPUTS:
%   x_tx                Transmitter positions [m]
%   x_rx                Receiver positions [m]
%   x_tgt               Target positions [m]
%   v_tx                Transmitter velocities [m/s]
%   v_rx                Receiver velocities [m/s]
%   v_tgt               Target velocities [m/s]
%   ref_idx             Matrix of tx/rx pairing indices (in the order
%                       they're used in C).  If ignored, then all pairwise
%                       measurements are used (nTx x nRx)
%
% OUTPUTS:
%   J               nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position
%
% Nicholas O'Donoughue
% 10 September 2021

% Parse inputs
[nDim,nTx] = size(x_tx);
[nDim2,nRx] = size(x_rx);
[nDim3,nTgt] = size(x_tgt);

assert(nDim == nDim2 && nDim2 == nDim3, ...
       'Input variables must match along first dimension.');

if nargin < 7 || ~exist('ref_idx','var') || isempty(ref_idx)
    ref_idx = [];
end

epsilon = 1e-15; % minimum range; to avoid divide by zero errors

% Compute Tx-Tgt Ranges
dxt = reshape(x_tgt,nDim,1,nTgt) - reshape(x_tx,nDim,nTx);
Rt = max(sqrt(sum(abs(dxt).^2,1)),epsilon); % 1 x nTx x nTgt
Jt = dxt./Rt; % Tx-tgt component of jacobian

% Compute Rx-Tgt Ranges
dxr = reshape(x_tgt,nDim,1,nTgt) - reshape(x_rx,nDim,nRx);
Rr = max(sqrt(sum(abs(dxr).^2,1)),epsilon); % 1 x nRx x nTgt
Jr = dxr./Rr; % Rx-Tgt component of jacobian

% Parse Tx/Rx pairing indices
if isempty(ref_idx) || strcmpi(ref_idx,'full')==0
    % Nothing (or full set of pairs) specified, do all pairs
    J = reshape(Jt,nDim,nTx,1,nTgt) + reshape(Jr,nDim,1,nRx,nTgt);
        % nDim x nTx x nRx x nTgt
        
    J = reshape(J,nDim, nTx*nRx, nTgt);
elseif isscalar(ref_idx)
    % Scalar specified, just keep that tx-rx pair
    J = Jt(:,ref_idx,:) + Jr(:,ref_idx,:);
else
    % Matrix; it must have two rows.  Each column is a tx/rx pair setting
    % transmitters indices on the first row, receiver on the second
    assert(size(ref_idx,1)==2,'Error parsing reference index; incorrect size.')
    assert(ismatrix(ref_idx),'Error parsing reference index; incorrect shape.');
    
    tx_idx = ref_idx(1,:);
    rx_idx = ref_idx(2,:);
    
    J = Jt(:,tx_idx,:) + Jr(:,rx_idx,:);
end
