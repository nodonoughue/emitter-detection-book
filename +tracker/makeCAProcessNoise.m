function Q = makeCAProcessNoise(maxG, num_states, heading, pitch, T, type)
% Q = makeCAProcessNoise(maxG, heading, pitch, T, type)
%
% Define the 9x9 process noise matrix for a 3D problem  with pos/vel/accel
% state space, with constant acceleration kinematic models, using one of
% three assumptions.
%
% type='lateral-only'
%       Air vehicle will divert laterally (x/y) with up to  maxG G's of 
%       lateral divert.  No thrusting changes.
%
% type='mostly-lateral'
%       Same as lateral only, except that 10% of divert process noise is
%       vertical, and 90% is lateral
%
% type='3d-divert'
%       Same as lateral only case, except the divert is allowed to
%       be in any dimension (x/y/z).
%
% INPUTS:
%   maxG        Maximum divert force (in G's)
%   heading     Heading (degrees CCW from +x axis)
%   pitch       Degrees above or below horizontal of target trajectory
%   T           Update rate (seconds)
%   type        String indicating the type of process noise to generate
%
% OUTPUTS:
%   Q           9 x 9 covariance matrix for process noise
%
% Nicholas O'Donoughue
% 11 Nov 2021


switch lower(type)
    case 'lateral-only'
        % No vertical acceleration, use heading along for lateral divert

        % 9x9 Process Noise Matrix
        qq_xyz = (9.81*maxG)^2/3 * [cosd(heading)^2 -cosd(heading)*sind(heading) 0;
                                -cosd(heading)*sind(heading) sind(heading)^2 0;
                                0 0 0];
    case '3d-divert'
        % No thrusting, only divert (in 3D)
        % Good for maneuvering aircraft under constant velocity
        u = [cosd(heading)*cosd(pitch), sind(heading)*cosd(pitch), sind(pitch)];
        P = eye(3) - u(:)*u(:)';
        qq_xyz = (9.81*maxG)^2/3 * P;
%        qq_xyz = (9.81*maxG)^2/3 * [cosd(heading)^2*cosd(pitch)^2, -cosd(heading)*sind(heading)*cosd(pitch)^2, cosd(heading)*sind(pitch)*cosd(pitch);
%                                -cosd(heading)*sind(heading)*cosd(pitch)^2, sind(heading)^2*cosd(pitch)^2, sind(heading)*cosd(pitch)*sind(pitch);
%                                cosd(heading)*cosd(pitch)*sind(pitch), sind(heading)*cosd(pitch)*sind(pitch), sind(pitch)^2];
    case 'mostly-lateral'
        % Allow 10% of divert in vertical (z) dimension, 90% in lateral
        % 9x9 Process Noise Matrix
        qq_xyz = (9.81*maxG)^2/3 * [cosd(heading)^2 -cosd(heading)*sind(heading) 0;
                                -cosd(heading)*sind(heading) sind(heading)^2 0;
                                0 0 0];
                            warning('mostly lateral case not properly defined');
end

% Pos/Vel/Accel distribution
j = [.5*T^2;T;1]; 
j = j(1:num_states); % keep only the first num_states states
J = j*j';

% Distribute x/y/z process noise to position, accel, and vel states
Q = kron(J,qq_xyz);