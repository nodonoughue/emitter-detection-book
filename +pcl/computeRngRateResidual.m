function err = computeRngRateResidual(x_tx,x_rx,x_tgt,v_tx,v_rx,rr_msrd)
% Compute the residual for a bistatic target at xtgt with unknown velocity,
% measured by a set of N transmitters at xtx with velocity vtx, and M 
% receivers at xrx with velocity vrx, given NxM range rate measurements 
% rr_msrd.
%
% Nicholas O'Donoughue
% 13 March 2020

% Parse inputs
[nDim,N] = size(x_tx);
[~,M] = size(x_rx);
[~,nTgt] = size(x_tgt);

% Compute distance from each source position to each sensor
d_tx = reshape(x_tgt,nDim,1,1,nTgt) - reshape(x_tx,nDim,N);
d_rx = reshape(x_tgt,nDim,1,1,nTgt) - reshape(x_rx,nDim,1,M);

Rtx = sqrt(sum(abs(d_tx).^2,1)); % 1 x nTx x nRx x nTgt
Rrx = sqrt(sum(abs(d_rx).^2,1)); % 1 x nTx x nRx x nTgt

% Reshape transmit/receive velocities
v_tx = reshape(v_tx,nDim,N);
v_rx = reshape(v_rx,nDim,1,M);

% Set up linear system of equations
%   y = NM x 1 vector of range rate measurements
%   m = NM x nDim matrix of bistatic range LOS vectors
%   b = NM x 1 vector of transmitter/receiver contributions to range rate

b = reshape(sum((d_tx.*v_tx./Rtx)+(d_rx.*v_rx./Rrx),1),N*M,nTgt); % NM x nTgt
m = reshape(-(d_tx./Rtx+d_rx./Rrx),nDim,N*M)'; % NM x nDim
y = reshape(rr_msrd,N*M,[]); % NM x nTgt

v_est = m \ (y-b); % nDim x nTgt solutions

% Compute residual
[~,rr] = pcl.measurement(x_tx,x_rx,x_tgt,v_tx,v_rx,v_est); % N x M x nTgt

err = reshape(rr,N*M,[])- y; % NM x nTgt