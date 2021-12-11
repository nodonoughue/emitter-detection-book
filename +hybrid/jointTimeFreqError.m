function C_r_rdot = jointTimeFreqError(a, sigma1, sigma2, delta_t, delta_f, f0, ts)
%


% Parse inputs
c = utils.constants.c;
kappa = delta_t/ts;
nu = 2*pi*delta_f*ts;

% Build the components for the FIM
m = -num_samples/2:(num_samples/2-1); % 13.13
Dv = diag(exp(-1i*nu*m)); % 13.12
Dk = diag(exp(-1i*2*pi*kappa/num_samples*m)); % 13.11
W = (1/sqrt(num_samples))*exp(-1i*2*pi/num_samples*(m*m')); % 13.10
M = diag(m);  % 13.18
Q = Dv*W'*Dk*W; % 13.19
overline_s = (2*pi/num_samples)*W'*M*W*s; % 13.16
s_tilde = Q*s; % 13.17

% Fisher Information Matrix
% Equation 13.15
const = 2/(a^2*sigma1^2+sigma2^2);
s_os = s'*overline_s;
stMst = s_tilde'*M*s_tilde;
Eos = overline_s'*overline_s;
osQMst = overline_s'*Q'*M*s_tilde;
stMMst = s_tilde'*M*M*s_tilde;

F = const * [Es, -s_os, stMst;
             -s_os, Eos, -real(osQMst);
             stMst, -real(osQMst), stMMst];

Finv = pinv(F);

% Build the covariance matrix
% Equation 13.21
H = [0, c*ts, 0;
     0, 0, c/(2*pi*f0*ts)];
C_r_rdot = H * Finv * H';