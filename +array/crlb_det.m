function C = crlb_det(Cs,sn2,psi_vec,M,v,v_dot)
% C = crlb_det(Cs,sn2,psi_vec,M,v,v_dot)
%
% Computes the CRLB for array-based DOA, according to
% section 9.4.
%
% INPUTS:
%
%   Cs              Source signal covariance matrix
%   sn2             Noise power
%   psi_vec         Array of D steering angles (in radians) for each source
%   M               Number of temporal snapshots taken
%   v               Function handle to steering vector v(psi)
%   v_dot           Function handle to steering vector gradient
%                   dv(psi)/dpsi
%
% OUTPUTS:
%   C               CRLB matrix (DxD) for angle of arrival estimation of
%                   each source, in radians--multiply by (180/pi)^2 to
%                   convert to degrees.
%
% Nicholas O'Donoughue
% 1 July 2019

% Build H Matrix
V = v(psi_vec);                     % N x D
D = v_dot(psi_vec);                 % N x D

Pv = (V / (V'*V))*V';               % N x N
Pv_ortho = eye(size(Pv)) - Pv;      % N x N
H = D'*Pv_ortho*D;                  % D x D

% Build the spectral matrix from linear SNR
Xi = Cs/sn2;                        % D x D

% CRLB
C = (1/(2*M)) * pinv(real(Xi.*H.'));        % D x D