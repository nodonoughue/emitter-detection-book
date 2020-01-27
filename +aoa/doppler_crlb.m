function crlb = doppler_crlb(snr,M,A,ts,f,R,fr,psi)
% crlb = doppler_crlb(snr,M,A,ts,f,R,fr,psi)
%
% Compute the lower bound on unbiased estimation error for a Doppler-based
% direction finding receiver.
%
% Inputs:
%
%   snr         Signal to Noise ratio [dB]
%   M           Number of samples
%   A           Received signal amplitude
%   ts          Sampling period [s]
%   f           Carrier frequency [Hz]
%   R           Doppler rotation radius [m]
%   fr          Doppler rotation rate [rotations/sec]
%   psi         True angle of arrival [radians]
%
% Outputs:
%
%   crlb        Lower bound on the Mean Squared Error of an unbiased
%               estimation of psi (radians)
%
% Nicholas O'Donoughue
% 1 July 2019


if numel(snr)>1
    crlb = reshape(arrayfun(@(x) aoa.doppler_crlb(x,M,A,ts,f,R,fr,psi),snr),size(snr));
    return;
end

% Build Jacobian of reference signal
snr_lin = 10.^(snr/10);
jr = 2*M.*snr_lin*[1./A^2 0 0 0;
                  0 1 ts.*(M-1)/2 0;
                  0 ts.*(M-1)/2 ts.^2.*(M-1).*(2*M-1)/6 0;
                  0 0 0 0];

% Build signal vectors and gradients
c = 3e8;
sx_dot = (ts*(0:M-1) + (R/c)*cos(2*pi*fr*ts*(0:M-1)-psi));
sx_ddot = 2*pi*f*(R/c)*sin(2*pi*fr*ts*(0:M-1)-psi);

% Build Jacobian of test (doppler) signal
Es_dot = sum(sx_dot);
Es_ddot = sum(sx_ddot);
Es_dot_dot = sum(sx_dot.^2);
Es_ddot_ddot = sum(sx_ddot.^2);
Es_dot_ddot = sum(sx_dot.*sx_ddot);

jx = 2*M*snr_lin*[1/A^2 0 0 0;
                  0 1 Es_dot/M Es_ddot/M;
                  0 Es_dot/M Es_dot_dot/M Es_dot_ddot/M;
                  0 Es_ddot/M Es_dot_ddot/M Es_ddot_ddot/M];

% Construct full jacobian and invert
J_full = jr + jx;
C_full = pinv(J_full);
crlb = C_full(end,end); % output in radians

