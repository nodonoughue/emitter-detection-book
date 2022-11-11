function J = grad_b(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source,tdoa_ref_idx,fdoa_ref_idx, do2DAoA, alpha_aoa, alpha_tdoa, alpha_fdoa)
% J = grad_b(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source,tdoa_ref_idx,...
%            fdoa_ref_idx, do2DAoA, alpha_aoa, alpha_tdoa, alpha_fdoa)
%
% Returns the gradient of hybrid measurements, with sensor uncertainties,
% with respect to sensor position and velocity.
%
%
% INPUTS:
%   x_aoa           AOA sensor positions
%   x_tdoa          TDOA sensor positions
%   x_fdoa          FDOA sensor positions
%   v_fdoa          FDOA sensor velocities
%   x_source        Candidate source positions
%   tdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for TDOA measurements
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA measurements
%   do2DAOA
%   alpha_aoa
%   alpha_tdoa
%   alpha_fdoa
%
% OUTPUTS:
%   J               nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position
%
% Nicholas O'Donoughue
% 31 October 2021

if nargin < 6 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

if nargin < 7 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

if nargin < 8 || ~exist('do2DAoA','var')
    do2DAoA = [];
end

if nargin < 9 || ~exist('alpha_aoa','var')
    alpha_aoa = [];
end

if nargin < 10 || ~exist('alpha_tdoa','var')
    alpha_tdoa = [];
end

if nargin < 11 || ~exist('alpha_fdoa','var')
    alpha_fdoa = [];
end

n_aoa = size(x_aoa,2);
n_tdoa = size(x_tdoa,2);
n_fdoa = size(x_fdoa,2);
n_source = size(x_source,2);

% Compute Jacobian for AOA measurements
if n_aoa>0
    J_aoa = triang.grad_b(x_aoa, x_source, do2DAoA, alpha_aoa);
else
    J_aoa = [];
end

% Compute Jacobian for TDOA measurements
if n_tdoa>0
    J_tdoa= tdoa.grad_b(x_tdoa, x_source, tdoa_ref_idx, alpha_tdoa);
else
    J_tdoa = [];
end

% Compute Jacobian for FDOA measurements
if n_fdoa>0
    J_fdoa= fdoa.grad_b(x_fdoa, v_fdoa, x_source, fdoa_ref_idx, alpha_fdoa);
else
    J_fdoa = [];
end

% Combine component Jacobians
if n_source <= 1
    % only one source pos, use blkdiag
    J = blkdiag(J_aoa,J_tdoa,J_fdoa);
else
    % build the jacobian by expanding along the first and third dimensions
    [n_beta_a,m_aoa,~] = size(J_aoa);
    [n_beta_t,m_tdoa,~] = size(J_tdoa);
    [n_beta_f,m_fdoa,~] = size(J_fdoa);

    J_1 = cat(2, J_aoa, zeros(n_beta_a,m_tdoa+m_fdoa,n_source));
    J_2 = cat(2, zeros(n_beta_t,m_aoa,n_source), J_tdoa, zeros(n_beta_t,m_fdoa,n_source));
    J_3 = cat(2,zeros(n_beta_f,m_aoa+m_tdoa,n_source),J_fdoa);
    J = cat(1,J_1,J_2,J_3);
end