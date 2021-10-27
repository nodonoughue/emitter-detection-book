function J = jacobian(x_tx, x_rx, x_tgt, v_tx, v_rx, v_tgt, ref_idx, angle_dims)
% J = jacobian(x_tx, x_rx, x_source, v_tx, v_rx, v_source, ref_idx, angle_dims)
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
%   angle_dims       Flag determining how many angle of arrival dimensions 
%                    to report (0 =  no angle data, 1 = azimuth, 2 =
%                    azimuth & elevation)
%
% OUTPUTS:
%   J               nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position
%
% Nicholas O'Donoughue
% 10 September 2021

%% Parse inputs
[nDim,nTx] = size(x_tx);
[nDim2,nRx] = size(x_rx);
[nDim3,nTgt] = size(x_tgt);

assert(nDim == nDim2 && nDim2 == nDim3, ...
       'Input variables must match along first dimension.');

if nargin < 8 || ~exist('angle_dims','var') || isempty(angle_dims)
    angle_dims = 0;
end
assert(angle_dims >=0 && angle_dims <= 2,'Error parsing angle dimensions command, must be between 0 and 2.');

if nargin < 7 || ~exist('ref_idx','var') || isempty(ref_idx)
    ref_idx = [];
end

do_doppler = nargin >= 4 && (~isempty(v_tx) || ~isempty(v_rx) || ~isempty(v_tgt));
if do_doppler && isempty(v_tx)
    v_tx = zeros(size(x_tx));
end

if do_doppler && isempty(v_rx)
    v_rx = zeros(size(x_rx));
end

if do_doppler && isempty(v_tgt)
    v_tgt = zeros(size(x_tgt));
end

%% Range Jacobian
epsilon = 1e-15; % minimum range; to avoid divide by zero errors

% Compute Tx-Tgt Ranges
dxt = reshape(x_tgt,nDim,1,nTgt) - reshape(x_tx,nDim,nTx);
Rt = max(sqrt(sum(abs(dxt).^2,1)),epsilon); % 1 x nTx x nTgt
Jt_rng = dxt./Rt; % Tx-tgt component of jacobian

% Compute Rx-Tgt Ranges
dxr = reshape(x_tgt,nDim,1,nTgt) - reshape(x_rx,nDim,nRx);
Rr = max(sqrt(sum(abs(dxr).^2,1)),epsilon); % 1 x nRx x nTgt
Jr_rng = dxr./Rr; % Rx-Tgt component of jacobian

%% Range-Rate Jacobian
if do_doppler
    u_tx = Jt_rng; % The transmitter-tgt jacobian is u_m(x) for transmitter m
    u_rx = Jr_rng; % The receiver-tgt jacobian is u_n(x) for receiver n
    
    % Projection matrix
    Proj_ortho_tx = eye(nDim) - reshape(u_tx,nDim,1,nTx,nTgt).*reshape(u_tx,1,nDim,nTx,nTgt);
    Proj_ortho_rx = eye(nDim) - reshape(u_rx,nDim,1,nRx,nTgt).*reshape(u_rx,1,nDim,nRx,nTgt);
    
    % Scaled velocity vector
    uv_t = (reshape(v_tx,nDim,nTx) - reshape(v_tgt,nDim,1,nTgt)) ./ Rt; % nDim x nTx x nTgt
    uv_r = (reshape(v_rx,nDim,nRx) - reshape(v_tgt,nDim,1,nTgt)) ./ Rr; % nDim x nRx x nTgt
    Jt_velrng = squeeze(sum(Proj_ortho_tx.*reshape(uv_t,1,nDim,nTx,nTgt),2)); % nDim x nTx x nTgt
    Jr_velrng = squeeze(sum(Proj_ortho_rx.*reshape(uv_r,1,nDim,nRx,nTgt),2)); % nDim x nRx x nTgt
    
    % Jacobian of range-rate w.r.t velocity is equivalent to the negative 
    % of the jacobian of range w.r.t. position
    Jt_vel = -Jt_rng;
    Jr_vel = -Jr_rng;
end

%% Angle of Arrival Measurements
if angle_dims > 0 
    do2DAOA = (angle_dims==2);
    
    Jr_psi = triang.jacobian(x_rx, x_tgt, do2DAOA);
        % This is occurring for each Tx/Rx pair, but triang is just doing
        % the geometry for target to receiver, so let's repeat it once
        % for each transmitter to capture the fact that each Tx/Rx link
        % has a unique look at that DOA.

    Jr_psi = reshape(Jr_psi,nDim,1,nRx,angle_dims,nTgt);
else
    Jr_psi = [];
end

%% Parse Tx/Rx pairing indices
if isempty(ref_idx) || strcmpi(ref_idx,'full')==1
    % Nothing (or full set of pairs) specified, do all pairs
    Jrng = reshape(Jt_rng,nDim,nTx,1,nTgt) + reshape(Jr_rng,nDim,1,nRx,nTgt);
    if angle_dims > 0
        Jpsi = reshape(repmat(Jr_psi,1,nTx,1,1,1),nDim,nTx,nRx*angle_dims,nTgt);
    end
    if do_doppler
        Jrng = cat(1,Jrng,zeros(size(Jrng))); % the jacobian of range w.r.t. vel is zero
        if angle_dims > 0
            Jpsi = cat(1,Jpsi,zeros(size(Jpsi))); % the jacobian of DOA w.r.t. vel is zero
        end
        Jvelrng = reshape(Jt_velrng,nDim,nTx,1,nTgt) + reshape(Jr_velrng,nDim,1,nRx,nTgt);
        Jvel = reshape(Jt_vel,nDim,nTx,1,nTgt) + reshape(Jr_vel,nDim,1,nRx,nTgt);
            % Jrng is (2*nDim) x nTx x nRx x nTgt
            % Jpsi is (2*nDim) x nTx x (nRx*angleDims) x nTgt
            % Jvelrng and Jvel are nDim x nTx x nRx x nTgt
    end
    
    if do_doppler
        if angle_dims > 0
            J = cat(2,reshape(Jrng,2*nDim,nTx*nRx,nTgt), reshape(cat(1,Jvelrng,Jvel),2*nDim,nTx*nRx,nTgt), reshape(Jpsi,2*nDim,nTx*nRx*angle_dims,nTgt));
        else
            J = cat(2,reshape(Jrng,2*nDim,nTx*nRx,nTgt), reshape(cat(1,Jvelrng,Jvel),2*nDim,nTx*nRx,nTgt));
        end
            % (2*nDim) x ((2+angle_dims)*nTx*nRx) x nTgt
            % nMsmt = (2+angle_dims)*nTx*nRx
    else
        if angle_dims > 0
            J = cat(2,reshape(Jrng,nDim, nTx*nRx, nTgt), reshape(Jpsi,nDim, nTx*nRx*angle_dims, nTgt));
        else
            J = reshape(Jrng,nDim, nTx*nRx, nTgt);
        end
            % nDim x ((1+angle_dims)*nTx*nRx) x nTgt
            % nMsmt = (1+angle_dims)*nTx*nRx
    end
elseif isscalar(ref_idx)
    % Scalar specified, just keep that tx-rx pair
    Jrng = Jt_rng(:,ref_idx,:) + Jr_rng(:,ref_idx,:);
    if angle_dims > 0
        Jpsi = reshape(Jr_psi(:,1,ref_idx,:,:),nDim,angle_dims,nTgt);
    end
    if do_doppler
        Jrng = cat(1,Jrng,zeros(size(Jrng)));
        if angle_dims > 0
            Jpsi = cat(1,Jpsi,zeros(size(Jpsi)));
        end
        Jvelrng = Jt_velrng(:,ref_idx,:) + Jr_velrng(:,ref_idx,:);
        Jvel = Jt_vel(:,ref_idx,:) + Jr_vel(:,ref_idx,:);
        if angle_dims > 0
            J = cat(2,Jrng,cat(1,Jvelrng,Jvel),Jpsi);
        else
            J = cat(2,Jrng,cat(1,Jvelrng,Jvel));
        end
            % (2*nDim) x (2+angle_dims) x nTgt
            % nMsmt = 2 + angle_dims (single tx/rx pair)
    else
        if angle_dims > 0
            J = cat(2,Jrng,Jpsi);
                % nDim x (1 + angle_dims) x nTgt
                % nMsmt = 1 + angle_dims
        else
            J = Jrng;
                % nDim x 1 x nTgt
                % nMsmt = 1
        end
    end
else
    % Matrix; it must have two rows.  Each column is a tx/rx pair setting
    % transmitters indices on the first row, receiver on the second
    assert(size(ref_idx,1)==2,'Error parsing reference index; incorrect size.')
    assert(ismatrix(ref_idx),'Error parsing reference index; incorrect shape.');
    
    tx_idx = ref_idx(1,:);
    rx_idx = ref_idx(2,:);
    
    Jrng = Jt_rng(:,tx_idx,:) + Jr_rng(:,rx_idx,:);
    if angle_dims > 0
        Jpsi = reshape(Jr_psi(:,1,rx_idx,:,:),nDim,numel(rx_idx)*angle_dims,nTgt);
    end
    if do_doppler
        Jrng = cat(1,Jrng,zeros(size(Jrng)));
        if angle_dims > 0
            Jpsi = cat(1,Jpsi,zeros(size(Jpsi)));
        end
        Jvelrng = Jt_velrng(:,tx_idx,:) + Jr_velrng(:,rx_idx,:);
        Jvel = Jt_vel(:,tx_idx,:) + Jr_vel(:,rx_idx,:);

        if angle_dims > 0
            J = cat(2,Jrng,cat(1,Jvelrng,Jvel), Jpsi);
                % (2*nDim) x (2+angle_dims)*size(ref_idx,2) x nTgt
                % nMsmt = (2+angle_dims)*size(ref_idx,2)
        else
            J = cat(2,Jrng,cat(1,Jvelrng,Jvel));
                % (2*nDim) x 2*size(ref_idx,2) x nTgt
                % nMsmt = 2*size(ref_idx,2)
        end
    else
        if angle_dims > 0
            J = cat(2,Jrng,Jpsi);
                % nDim x (1+angle_dims)*size(ref_idx,2) x nTgt
                % nMsmt = (1+angle_dims)*size(ref_idx,2)
        else
            J = Jrng;
                % nDim x size(ref_idx,2) x nTgt
                % nMsmt = size(ref_idx,2)
        end
    end
end