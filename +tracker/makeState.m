function s = makeState(state_space, time, state, covar)
% makeState  Create a tracker state struct.
%
% s = makeState(state_space, time)
% s = makeState(state_space, time, state)
% s = makeState(state_space, time, state, covar)
%
% INPUTS
%   state_space  Struct from makeMotionModel,
%                with at minimum the fields:
%                  num_states  total length of the state vector
%                  num_dims    number of spatial dimensions
%                  has_vel     logical, true if velocity is a tracked state
%                  pos_idx     index vector for position elements
%                  vel_idx     index vector for velocity elements (if has_vel)
%   time         Scalar timestamp [s]
%   state        Column vector of length state_space.num_states.
%                Default: zeros(num_states, 1)
%   covar        Square covariance matrix (num_states x num_states).
%                Default: [] (no covariance)
%
% OUTPUTS
%   s   Struct with fields:
%         state_space  – the state_space struct (shared reference)
%         time         – scalar timestamp [s]
%         state        – (num_states x 1) state vector
%         covar        – (num_states x num_states) covariance matrix, or []
%
% Nicholas O'Donoughue
% June 2025

if nargin < 3 || isempty(state)
    state = zeros(state_space.num_states, 1);
end

if nargin < 4
    covar = [];
end

s = struct('state_space', state_space, ...
           'time',        time, ...
           'state',       state(:), ...
           'covar',       covar);
