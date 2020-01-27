function psi = doppler_df(r,x,ts,f,R,fr,psi_res,min_psi,max_psi)
% psi = doppler_df(r,x,ts,f,R,fr,psi_res,min_psi,max_psi)
%
% Compute the estimate of a signal's angle of arrival (in radians), given a
% reference signal, and Doppler signal (from an antenna rotating about the
% reference antenna).
%
% Inputs:
%
%   r           Reference signal, length M
%   x           Doppler test signal, length M
%   ts          Sampling period [s]
%   f           Carrier frequency [Hz]
%   R           Doppler rotation radius [m]
%   fr          Doppler rotation rate [rotations/sec]
%   psi_res     Desired output AoA resolution [radians]
%   min_psi     Minimum bound on valid region for psi [default = -pi]
%   max_psi     Maximum bound on valid region for psi [default = pi]
%
% Outputs:
%
%   psi         Estimated angle of arrival
%
% Nicholas O'Donoughue
% 1 July 2019

% Speed of light
c=3e8;

% Filter the reference signal out of the test signal
y = x.*conj(r);

% Compute (and unwrap) the complex phase angle
phi = unwrap(atan2(imag(y),real(y)));

% Generate test phase signal
t_vec = ts*(0:numel(phi)-1);
phi_0 = @(psi) (2*pi*f*R/c)*cos(2*pi*fr*t_vec(:) - psi(:)');

% Initialize DF search loop
this_psi_res = 1;                       % Start at 1 radian resolution
psi_vec = min_psi:this_psi_res:max_psi; % Set up search vector

% Loop until desired DF resolution achieved
while this_psi_res > psi_res
    % Compute error at each test point in search vector
    err = sum(abs(phi(:)-phi_0(psi_vec)).^2,1);
    
    % Find point with minimal error, use as center point for next
    % iteration of search
    [~,idx_opt] = min(err);
    psi = psi_vec(idx_opt);
    
    % Increase resolution (descrease search vector step size)
    % and interval; create search vector for next iteration
    this_psi_res = this_psi_res/10;
    min_psi = psi_vec(max(1,idx_opt-2));
    max_psi = psi_vec(min(numel(psi_vec),idx_opt+2));
    psi_vec = min_psi:this_psi_res:max_psi;
end