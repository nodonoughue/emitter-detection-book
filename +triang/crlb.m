function C = crlb(x_aoa,x0,C_psi)
% C = crlb(x_aoa,x0,C_psi)
%
% Computes the CRLB for triangulation given the angle of arrival covariance
% matrix C_psi, sensor positions xs, and source position x0
%
% INPUTS:
%   x_aoa       2 x N list of sensor positions in 2D
%   x0          2 x M list of source positions in 2D
%   C_psi       Array covariance matrix (N x N) for angle of arrival
%               measurements
%
% OUTPUTS:
%   C           2 x 2 x M array of CRLB on source position estimates for
%               each source location
%
% Nicholas O'Donoughue
% 1 July 2019

% Loop over sensor positions
M = size(x0,2);
C = zeros(2,2,M);
for idx_m = 1:M
    x_m = x0(:,idx_m);
    
    % Handle edge cases
    R = utils.rng(x_aoa,x_m);
    bad_sensor = R < 1e-6;
    
    if sum(~bad_sensor) < 2
        % Insufficient sensors to provide ranging
        C(:,:,idx_m) = Inf*ones(2,2);
        continue;
    end
    
    % Compute Jacobian and remove bad sensors
    J = triang.jacobian(x_aoa,x_m); % 2 x N
    J = J(:,~bad_sensor);
    
    % Fisher Information Matrix
    F = J/C_psi(~bad_sensor,~bad_sensor)*J'; % 2 x 2
    
    % Invert the Fisher Information Matrix to compute the CRLB
    C(:,:,idx_m) = pinv(F); % 2 x 2
end

