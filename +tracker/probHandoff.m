function [P_handoff, C_proj] = probHandoff( C_track, psi_max, rng_max, el_min, el_max, bng_min, bng_max )
%PROBHANDOFF Summary of this function goes here
%   Compute probability of seeker handoff, given a track error covariance
%   specified by C_track.
%
%
% INPUTS:
%   C_track     Ndim x Ndim x M set of track error covariances
%   psi_max     Max cone angle between true and reported target position,
%               to support seeker handoff
%   rng_max     Max range at which seeker acqusition can occur
%   el_min      Minimum elevation angle for seeker on approach
%   el_max      Maximum elevation angle for seeker on approach
%   bng_min     Minimum bearing angle (deg CCW from +x) for seeker on 
%               approach
%   bng_max     Maximum bearning angle (deg CCW from +x) for seeker on
%               approach
%
% OUTPUTS:
%   P_handoff   M x 1 vector of probabilities for seeker handoff, averaged
%               over the intervals bng_min:bng_max for bearing and
%               el_min:el_max for elevation.
%


%% Parse Inputs
assert(nargin >= 4, 'Not enough input arguments.');
[Nd1, Nd2, M] = size(C_track);
assert(Nd1==Nd2 && Nd1==3, 'First two dimensions of C_track must have 3 elements.');

assert(~isempty(el_min) || (nargin >= 5 && ~isempty(el_max)),...
        'At least one elevation angle for seeker approach must be defined.');
assert((nargin >= 6 && ~isempty(bng_min)) || (nargin >= 7 && ~isempty(bng_max)),...
        'At least one bearing angle for seeker approach must be defined.');

%% Process Cone Angle and Distance
max_radius = tand(psi_max)*rng_max;

%% Build Vectors
if isempty(el_max)
    el_vec = el_min;
elseif isempty(el_min)
    el_vec = el_max;
else
    el_vec = linspace(el_min,el_max,101);
end
N_el = numel(el_vec);

if nargin < 7 || isempty(bng_max)
    bng_vec = bng_min;
elseif isempty(bng_min)
    bng_vec = bng_max;
else
    bng_vec = linspace(bng_min,bng_max,101);
end
N_bng = numel(bng_vec);

%% Loop over Covariance matrices
P_handoff = zeros(M,1);
C_proj = zeros(3, 3, M);

for m=1:M
    this_C = squeeze(C_track(:,:,m));
    
    if any(isnan(this_C))
        warning('Track error covariance has NaN values, unable to compute Probability of Handoff.');
        P_handoff(m) = NaN;
        continue;
    end
    % Loop over az/el settings
    this_P = zeros(N_el,N_bng);
    for i_bng = 1:N_bng
        this_bng = bng_vec(i_bng);
        cos_bng_i = cosd(this_bng);
        sin_bng_i = sind(this_bng);
        
        for i_el = 1:N_el
            this_el = el_vec(i_el);
            cos_el_i = cosd(this_el);
            sin_el_i = sind(this_el);
            
             % Computer seeker position at acquisition
             x_seeker = rng_max * [cos_bng_i*cos_el_i;
                                   sin_bng_i*cos_el_i;
                                   sin_el_i];
            
             % Project Error
             this_C_proj = utils.projectError(this_C,[0;0;0],x_seeker);
             
             % Ellipse CDF
             this_P(i_el,i_bng) = utils.ellipseCDF(this_C_proj, max_radius);
             
             if isnan(this_P(i_el, i_bng))
                 warning('NaN probability encountered');
             end
            if i_bng == 1 && i_el==1
                C_proj(:,:,m) = this_C_proj;
            end
        end
    end
    
    P_handoff(m) = sum(this_P(:))/(N_el*N_bng);
end

