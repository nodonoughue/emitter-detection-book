function x = chanHoSoln(x0,rho,C,ref_idx)
% x = chanHoSoln(x0,rho,C,ref_idx)
%
% Computes the Chan-Ho solution for TDOA processing.
%
% Inputs:
%   
%   x0              Sensor positions [m]
%   rho             Range-Difference Measurements [m]
%   C               Error covariance matrix for Range-Difference
%                   measurements [m^2]
%   ref_idx         Reference sensor index
%
% Outputs:
%   x               Estimated source position [m]
%
% Ref:  Y.-T. Chan and K. Ho, “A simple and efficient estimator for 
%       hyperbolic location,” IEEE Transactions on signal processing, 
%       vol. 42, no. 8, pp. 1905–1915, 1994.
%
% Nicholas O'Donoughue
% 1 July 2019

[num_dims,num_sensors] = size(x0);
if nargin < 4 || isempty(ref_idx)
    ref_idx = num_sensors;
end

%% Stage 0: Parse Covariance Matrix
% The Chan-Ho solution implicitly uses the final sensor as a common
% reference
test_idx = setdiff(1:num_sensors,ref_idx);
C = utils.resampleCovMtx(C, test_idx, ref_idx);

%% Stage 1: Initial Position Estimate

% Compute system matrix overline(A) according to 13.23

% Compute shifted measurement vector overline(y) according to 13.24
R = sqrt(sum(abs(x0).^2,1))';
y1 = (rho.^2-R(test_idx).^2+R(ref_idx)^2);
G1 = -2*[(x0(:,test_idx) - x0(:,ref_idx))',rho];

% Compute initial position estimate overline(theta) according to 13.25
B = eye(num_sensors-1);
W1 = (B*C*B');
th1 = ((G1'/W1)*G1)\(G1'/W1)*y1;

% Refine sensor estimate
for j=1:3
    ri_hat = sqrt(sum((x0-th1(1:end-1)).^2,1));
    B = 2*diag(ri_hat(1:end-1));
    W1 = (B*C*B');
    th1 = ((G1'/W1)*G1)\(G1'/W1)*y1;
end

th1p = th1 - [x0(:,ref_idx);0];


%% Stage 2: Account for Parameter Dependence

y2 = th1.^2;
G2 = cat(1,eye(num_dims),ones(1,num_dims));

B2 = 2*diag(th1p);

% Compute final parameter estimate overline(theta)' according to 13.32
W2 = B2'\(G1'*W1*G1)/B2;
th2 = (G2'*W2*G2)\G2'*W2*y2;

% Compute position estimate overline(x)' according to 13.33
% x_prime1 = x0(:,end)+sqrt(th_prime);
% x_prime2 = x0(:,end)-sqrt(th_prime);
% 
% offset1 = norm(x_prime1-th(1:end-1));
% offset2 = norm(x_prime2-th(1:end-1));
% 
% if offset1 <= offset2
%     x = x_prime1;
% else
%     x = x_prime2;
% end
x = sign(diag(th1(1:end-1)))*sqrt(th2)+x0(:,ref_idx);