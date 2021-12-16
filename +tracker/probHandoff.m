function P_handoff = probHandoff( C_track, psi_max, rng_max, el_min, el_max, bng_min, bng_max )
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
max_radius = tan(psi_max)*rng_max;

%% Build Vectors
if isempty(el_max)
    el_vec = el_min;
elseif isempty(el_min)
    el_vec = el_max;
else
    el_vec = linspace(el_min,el_max,101);
end
N_el = numel(el_vec);

if isempty(bng_max)
    bng_vec = bng_min;
elseif isempty(bng_min)
    bng_vec = bng_max;
else
    bng_vec = linspace(bng_min,bng_max,101);
end
N_bng = numel(bng_vec);

%% Loop over Covariance matrices
P_handoff = zeros(M,1);
for m=1:M
    this_C = squeeze(C_track(:,:,m));
    
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
             C_proj = utils.projectError(this_C,[0;0;0],x_seeker);
             
             % Compute principal eigenvalues
             [~,lam] = eig(C_proj);
             lamSort = sort(lam,'descend');
             sigma2 = lamSort(1:2); % Variance of the two principal components
             
             % Compute probability that a random draw with these variances
             % will have radius < max_radius.
             %
             % The radius of error follows a Hoyt distribution (special
             % case of Nakagami distribution)
             % https://reference.wolfram.com/language/ref/HoytDistribution.html
             %
             % Formulated differently, the square of the radius of error is
             % the weighted sum of two chi-squared random variables with
             % one degree of freedom (k=1), and weights given by sigma.
             % See https://arxiv.org/pdf/1208.2691.pdf, Corollary 1
             %
             % In the meantime, let us approximate:
             %   1. Find the volume of an ellipse such that the radius of
             %      the larger axis = max_radius.  This will be a lower
             %      bound on P{Handoff}
             %   2. Find the volume of an ellipse such that the radius of
             %      the smaller axis = max_radius. This will be a lower
             %      bound on P{Handoff}
             %   3. Average the two results
             
             gamma = max_radius^2 ./ sigma2; % Scale factor to make each component's variation = the max acceptable error
             prob_handoff_bounds = utils.computeRMSEConfInterval(gamma);
             this_P(i_el,i_bng) = mean(prob_handoff_bounds);
        end
    end
    
    P_handoff(m) = sum(this_P(:))/(N_el*N_bng);
end

