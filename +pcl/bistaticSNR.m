function snr = bistaticSNR(x_tx,x_rx,x_tgt,erp,bt,f0,gr,rcs,bw,nf,loss,ref_idx)
% snr = bistaticSNR(x_tx,x_rx,x_tgt,erp,bt,f0,gr,rcs,bw,l,ref_idx)
%
% Computes the bi-static SNR, according to equation 2.20 from
% Malanowski, Bistatic Radar Signal Processing, Artech House, 2019
%
% This is the *output* SNR (after correlation processing).
%
% The required size of all non-singleton dimensions is provided below.
%
% x_tx          Transmitter positions (ndim x ntx)
% x_rx          Receiver positions (ndim x nrx)
% x_tgt         Target positions (ndim x ntgt)
% erp           Transmitter ERP [dBW] (ntx x 1)
% bt            Time-Bandwidth Product [unitless] (ntx x 1)
% f0            Carrier Frequency [Hz] (ntx x 1)
% gr            Receive antenna gain [dBi] (ndim x nrx x ntgt)
% rcs           Bistatic RCS [dBsm] (ntx x nrx x ntgt)
% bw            Receiver bandwidth [Hz] (ntx x nrx)
% nf            Noise Figure [dB] (nrx x 1)
% loss          Receiver loss [dB] (nrx x 1)
% ref_idx       Scalar, or 2xnMsmt vector of tx/rx indices to use
%
% OUTPUT:
%   snr         Signal-to-Noise Ratio [dB] (1 x nMsmt)
%
% Nicholas O'Donoughue
% 26 October 2021


ntx = size(x_tx,2);
nrx = size(x_rx,2);
ntgt = size(x_tgt,2);

% Compute Ranges
Rtx = utils.rng(x_tx,x_tgt); % ntx x ntgt
Rrx = utils.rng(x_rx,x_tgt); % nrx x ntgt

% Wavelength
lambda = utils.constants.c/f0;

% Tx FOM
% ERP * bt * (lam^2) / (4*pi)^3 Rtx^2 
tx_fom = erp + 10*log10(bt) + 20*log10(lambda) - 30*log10(4*pi) - 20*log10(Rtx);
    % size should be ntx x ntgt

% Rx FOM
% Gr / Rrx^2 * N * L
N = 10*log10(utils.constants.kT*bw) + nf;
rx_fom = gr - 20*log10(Rrx) - N - loss;

%% Parse ref_idx
if isempty(ref_idx) || strcmpi(ref_idx,'full')==1
    % Nothing specified, or the command 'full', do all tx/rx sensor
    % pairings
    snr = reshape(tx_fom,ntx,1,ntgt) + reshape(rx_fom,1,nrx,ntgt) + rcs;

    snr = reshape(snr,ntx*nrx,ntgt);
elseif isscalar(ref_idx)
    % Scalar specified, keep that pair

    % Parse RCS
    if size(rcs,1) > 1 % indexed across transmitter
        tx_idx = ref_idx;
    else
        tx_idx = 1;
    end
    if size(rcs,2) > 1 % indexed across receiver
        rx_idx = ref_idx;
    else
        rx_idx = 1;
    end
    this_rcs = rcs(tx_idx,rx_idx,:);

    snr = tx_fom(ref_idx,:) + rx_fom(ref_idx,:) + this_rcs;

else
    % Matrix: it must have two rows, one for tx ID and one for rx ID
    assert(size(ref_idx,1)==2,'Error parsing reference index; incorrect size.')
    assert(ismatrix(ref_idx),'Error parsing reference index; incorrect shape.');
    
    tx_idx = ref_idx(1,:);
    rx_idx = ref_idx(2,:);

    % Parse RCS
    if size(rcs,1)>1
        rcs_tx_idx = tx_idx;
    else
        rcs_tx_idx = 1;
    end
    if size(rcs,1)>1
        rcs_rx_idx = rx_idx;
    else
        rcs_rx_idx = 1;
    end

    % Construct SNR
    snr = tx_fom(tx_idx,:) + rx_fom(rx_idx,:) + rcs(rcs_tx_idx,rcs_rx_idx,:);
            % size(ref_idx,2) x ntgt
end