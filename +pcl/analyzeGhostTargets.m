function ghost_tgts = analyzeGhostTargets(xtx,xrx,xtgt,vtgt,rng_std_dev,rr_std_dev,doa_std_dev)
% ghost_tgts = analyzeGhostTargets(xtx,xrx,xtgt,vtgt,rng_std_dev,rr_std_dev,doa_std_dev)
%
% Perform PCL geolocation for a set of transmitters Xtx and receivers Xrx,
% with bistatic measurements msmt_struct, using the specified solver.
%
% Performs pair-wise matching of bistatic measurements, removing infeasible
% results.
%
% INPUTS:
%   xtx         Transmitter positions
%   xrx         Reciever positions
%   xtgt        Target positions
%   vtgt        Target velocities
%   rng_std_dev Bistatic range measurement error [default = 1000 m]
%   rr_std_dev  Bistatic range rate measurement error [default = .5 m/s]
%   doa_std_dev Receiver AOA estimation error [default = 360 deg]
%
% OUTPUTS:
%   tgt_struct  Array of detected target positions
%
%
% Nicholas O'Donoughue
% 25 March 2020

% Parse inputs
if nargin < 5 || isempty(rng_std_dev)
    rng_std_dev = 1e3; % use 1 km by default
end

if nargin < 6 || isempty(rr_std_dev)
    rr_std_dev = .5; % use .5 m/s by default
end

if nargin < 7 || isempty(doa_std_dev)
    doa_std_dev = Inf; % use INF by default (no DOA)
    do_doa =false;
else
    do_doa = isfinite(doa_std_dev);
end

% Initialize ghost target struct
n_ghost = 100;
idx_ghost = 1;
blank_ghost_struct = struct('tx_ids',[],'rx_id',[],'tgt_ids',[],'pos',[],'rng_err',[],'rr_err',[]);
ghost_tgts(n_ghost) = blank_ghost_struct;

% Compute measurements for receiver at the origin
[rng_bi, rr_bi] = pcl.measurement(xtx,xrx,xtgt,zeros(size(xtx)),zeros(size(xrx)),vtgt);
[psi,psi_elev] = triang.measurement(xrx,xtgt);

ntx = size(xtx,2);
nrx = size(xrx,2);
ntgt = size(xtgt,2);

% Look for all intersection points for 3 or more spheroids
for rx = 1:nrx
    this_x_rx = xrx(:,rx);
    for tx1 = 1:ntx
        for tx2 = 1:ntx
            if tx2==tx1, continue, end
            
            for tx3 = 1:ntx
                if tx3==tx2 || tx3==tx1, continue, end
                
                for tx4 = 1:ntx
                    if tx4==tx1 || tx4==tx2 || tx4==tx3, continue, end
                    tx_ids = [tx1,tx2,tx3,tx4];
                    this_x_tx = xtx(:,tx_ids);
                    this_x_tx_centered = this_x_tx - this_x_rx;
                    
                    % Look for all possible intersections
                    for tgt1 = 1:ntgt
                        for tgt2 = 1:ntgt
                            for tgt3 = 1:ntgt
                                for tgt4=1:ntgt
                                    
                                    if tgt1==tgt2 && tgt2==tgt3 && tgt3==tgt4, continue, end
                                    
                                    tgt_ids = [tgt1,tgt2,tgt3,tgt4];
                                    
                                    % Grab the bistatic range and doppler
                                    % measurements
                                    tmp = squeeze(rng_bi(:,rx,:));
                                    this_rng_bi = tmp(sub2ind([ntx,ntgt], tx_ids, tgt_ids))';
                                    tmp = squeeze(rr_bi(:,rx,:));
                                    this_rr_bi =  tmp(sub2ind([ntx,ntgt], tx_ids, tgt_ids))';
                                    
                                    this_psi = psi(rx,tgt_ids)';
                                    this_psi_elev = psi_elev(rx,tgt_ids)';
                                    
                                    % Compute solution
                                    x_ghost = pcl.sxSoln(this_x_tx_centered,this_rng_bi) + this_x_rx;
                                    
                                    % Skip if any of the three bistatic range msmts
                                    % are >2 std deviations from true
                                    r_ghost = pcl.measurement(this_x_tx,this_x_rx,x_ghost);
                                    err_rng_bi = abs(this_rng_bi - r_ghost);
                                    
                                    if any(err_rng_bi > 2*rng_std_dev)
                                        continue;
                                    end
                                    
                                    % Skip if any of the three bistatic range
                                    % rate measurements are >2 std deviations from
                                    % true
                                    err_rr_bi = pcl.computeRngRateResidual(this_x_tx,this_x_rx,x_ghost,zeros(size(this_x_tx)),zeros(3,1),this_rr_bi);
                                    
                                    if any(err_rr_bi > 2*rr_std_dev)
                                        continue
                                    end
                                    
                                    % Check angle of arrival
                                    if do_doa
                                        [psi_ghost, psi_elev_ghost] = triang.measurement(this_x_rx,x_ghost);
                                        
                                        psi_err = abs(this_psi - psi_ghost);
                                        psi_elev_err = abs(this_psi_elev - psi_elev_ghost);
                                        
                                        if any(psi_err > 2*doa_std_dev) || any(psi_elev_err > 2*doa_std_dev)
                                            continue
                                        end
                                    end
                                    
                                    % This is a ghost target, add it to the set
                                    if idx_ghost > n_ghost
                                        % Ghost target struct is full, expand it
                                        clear tmp
                                        tmp(2*n_ghost) = blank_ghost_struct;
                                        tmp(1:n_ghost) = ghost_tgts;
                                        ghost_tgts = tmp;
                                        clear tmp;
                                        
                                        n_ghost = numel(ghost_tgts);
                                    end
                                    
                                    ghost_tgts(idx_ghost) = struct('tx_ids',tx_ids,...
                                        'rx_id',rx,...
                                        'tgt_ids',tgt_ids,...
                                        'pos',x_ghost,...
                                        'rng_err',err_rng_bi,...
                                        'rr_err',err_rr_bi);
                                    % Increment counter
                                    idx_ghost = idx_ghost + 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

ghost_tgts = ghost_tgts(1:idx_ghost-1);