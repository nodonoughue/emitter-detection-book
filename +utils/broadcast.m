function [common_size, error_msg, varargout] = broadcast(varargin)
% [common_size, error_msg, ...] = broadcast(A1,A2,...)
%
% Checks the size of a set of N-dimensional inputs (A1, A2, ...).  If the
% inputs have a broadcastable size (e.g. all non-scalar dimensions match),
% then the common size is computed and returned.
%
% If there is a mismatch, then the common_size = NaN, and the second output
% error_msg is populated with a string containing the error message.
%
% Optional: if the final input is a boolean that is True, then repmat is
% used to cast all of the inputs to the common size and returned after the
% error message.
%

error_msg = '';

% Parse the input for the flag to resize all inputs and return them
do_resize = false;
if islogical(varargin{end}) && varargin{end}
    do_resize = true;
    varargin = varargin(1:end-1);
end

% Loop through inputs and compute common size
common_size = size(varargin{1});
for idx_in = 2:numel(varargin)
    this_sz = size(varargin{idx_in});
    
    % Find common size
    n_dim = numel(common_size);
    n_dim_new = numel(this_sz);
    n_dim_common = min(n_dim,n_dim_new);
    cmn_idx = 1:n_dim_common;
    
    % Compare first n_dim_common elements
    bad_match = any( (common_size(cmn_idx)>1) & (this_sz(cmn_idx)>1) & (common_size(cmn_idx)~=this_sz(cmn_idx)));
    if bad_match
        common_size = NaN;
        error_msg = 'Non-scalar input dimensions do not match.';
        varargout = {NaN};
        return;
    end
    
    % Replace any scalar entries in sz
    replace = common_size(cmn_idx)==1 & this_sz(cmn_idx)>1;
    common_size(replace) = this_sz(replace);
    
    % Add any new dimensions to sz
    if n_dim_new > n_dim
        new_idx = n_dim+1:n_dim_new;
        common_size(new_idx) = this_sz(new_idx);
    end
end


if do_resize
	varargout = cell(size(varargin));
    for idx_in = 1:numel(varargin)
        this_sz = size(varargin{idx_in});
        mult = common_size./this_sz;
        varargout{idx_in} = repmat(varargin{idx_in},mult);
    end
else
    varargout = {NaN};
end
