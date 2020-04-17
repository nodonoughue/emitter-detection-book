function x = sxSoln(Xtx,rng_bistatic)
% x = sxSoln(Xtx,Xrx,rng_bistatic)
%
% Computes the Spherical Intersection solution for Passive Coherent
% Location for a single receiver at the origin.
%
% Inputs:
%   
%   Xtx             Transmitter positions [m], Ndim x Ntx
%   rng_bistatic    Bistatic range measurements [m], Ntx x 1
%
% Outputs:
%   x               Estimated source positions [m], Ndim x 1
%
% Ref:  Mateusz Malanowski, Signal Processing for Passive Bistatic Radar,
%       Artech House, 2019. Chapter 8
%
% Nicholas O'Donoughue
% 25 March 2020

%% Stage 1 - Build matrices

z = .5 *  (sum(Xtx.^2,1)' - rng_bistatic.^2); % equation 8.10
a = Xtx'\z;  % equation 8.17
b = Xtx'\rng_bistatic;  % equation 8.18

%% Step 2 - Estimate target range from the receiver
% Equation 8.21

ab = a'*b;
aa = a'*a;
bb = b'*b;

Rt1 = real(-ab + sqrt(ab^2-(bb-1)*aa))/(bb-1);
Rt2 = real(-ab - sqrt(ab^2-(bb-1)*aa))/(bb-1);

%% Step 3 - Estimate target position
% Substitute Rt1 and Rt2 into equation 8.12,
% take the solution with the higher altitude

if Rt1 <= 0
    % Infeasible solution
    xtgt1 = [0 0 0]';
else
    xtgt1 = a + b*Rt1;
end

if Rt2 <= 0
    % Infeasible solution
    xtgt2 = [0 0 0]';
else
    xtgt2 = a + b*Rt2;
end

if xtgt1(3) > xtgt2(3)
    % Rt1 yields a higher altitude estimate, select that one
    x = xtgt1;
else
    % Rt2 yields a higher altitude estimate, select that one
    x = xtgt2;
end