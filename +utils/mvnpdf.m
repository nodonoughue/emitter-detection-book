function p = mvnpdf(x, mu, sigma)
% p = mvncdf(x, mu, sigma)
%
% Reduced functionality calculation of the PDF of a multi-variate random
% normal distribution with mean vector mu and covariance matrix sigma.
% To be used if the user does not have the Statistics & Machine Learning 
% Toolbox installed.
%
%
% INPUTS:
%   x           The value of the vector x at which to compute the PDF; each
%               row is an independent observation
%   mu          Expectation of x; must match the size of the rows of x. If
%               blank, then expectation is assumed to be zero for all
%               dimensions of x.
%   sigma       Covariance matrix of x (if blank, an identity matrix is
%               assumed); must match the size of the rows of x.
%
% OUTPUTS:
%   p           PDF for a multi-variate normal distribution with
%               expectation mu and covariance matrix sigma, evaluated at 
%               the values indicated by the columns of x.
%
% Nicholas O'Donoughue
% 13 April 2022

% Check Inputs
assert(nargin>=1,'Not enough inputs.');
assert(ismatrix(x),'Invalid data; observations must along the rows of a 2D matrix.')

% Parse Dimensions
if iscolumn(x), x = x.'; end  % If it's a column vector, make it a row vector
[num_obs, num_dim] = size(x);

% Compute X Offset
if nargin < 2 || isempty(mu)
    x0 = x;
else
    assert(ndims(mu)<=2,'Invalid data; expectation mu has unexpected dimensions.');
    if iscolumn(mu), mu = mu.'; end % If mu is a column vector, make it a row vector
    [r, c] = size(mu);

    % Case 1 -- mu is a row or column vector
    if isrow(mu)
        x0 = x - mu;
    % Case 2 -- separate mu for each observation
    else
        assert(c == num_dim, 'Invalid data; expectation mu has unexpected dimensions.');
        assert(num_obs == 1 || r == 1 || num_obs == r, 'Invalid data; expectation mu has unexpected dimensions.');
        
        x0 = x - mu;
        num_obs = max(num_obs, r); % Make sure num_obs matches the number of rows in x0
    end
end

% Compute intermediate data products
%  x_Rinv -- half of the quadratic term (x-x0)'*R^{-1}*(x-x0), one row
%            per observation
%  log_sqrt_det -- the log of the sqrt of the determinant of the covariance
%                   matrix (one per covariance matrix)
if nargin < 3 || isempty(sigma)
    % No covariance matrix provided
    x_Rinv = x0;
    log_sqrt_determinant = 0;
else
    % Possible covariance matrix input sizes
    % 1. Just the diagonal is supplied
    %       [1, num_dim]
    %       [1, num_dim, num_obs]
    % 2. Full matrix supplied
    %       [num_dim, num_dim]
    %       [num_dim, num_dim, num_obs]
    sz = size(sigma);
    if numel(sz) > 2 % multiple dimensions
        num_sigma = sz(3);
        assert(num_sigma==num_obs, 'Bad covariance matrix.');
    else
        num_sigma = 1;
    end
    
    if sz(1)==sz(2) % Full covariance matrix, single observation
        sigma_is_diag = false;
        assert(sz(1)==sz(2),'Bad covariance matrix.');
    else % Just the diagonal, should be 1 x num_dim
        sigma_is_diag = true;
        assert(sz(2)==num_dim,'Bad covariance matrix.');
    end

    if sigma_is_diag
        if num_sigma > 1
            R = reshape(sqrt(sigma),num_dim,num_sigma)'; % num_obs x num_dim
        else
            R = sqrt(sigma);
        end

        x_Rinv = x0 ./ R;
        log_sqrt_determinant = sum(log(R),2);
    else
        if num_sigma == 1
            R = chol(sigma,'upper');
            x_Rinv = x0 / R;
            log_sqrt_determinant = sum(log(diag(R)));
        else
            x_Rinv = zeros(num_obs, num_dim);
            log_sqrt_determinant = zeros(num_obs, 1);
            for idx_obs = 1:num_obs
                thisR = chol(squeeze(sigma(:,:,idx_obs)),'upper');
                x_Rinv(idx_obs, :) = x0(idx_obs, :) / thisR;
                log_sqrt_determinant(idx_obs) = sum(log(diag(thisR)));
            end
        end
    end 
end


% Compute the quadratic
quadratic = sum(x_Rinv.^2,2); % take the square, sum across dimensions

% PDF
p = exp(-0.5*quadratic - log_sqrt_determinant - num_dim*log(2*pi)/2);