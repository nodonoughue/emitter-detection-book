function xi = squareLawMinSNR(PFA,PD,M)
% xi = squareLawMinSNR(PFA,PD,M)
%
% Compute the required SNR to achieve the desired probability of detection,
% given the maximum acceptable probability of false alarm, and the number
% of complex samples M.
%
% The returned SNR is the ratio of signal power to complex noise power.
%
% Inputs:
%
%   PFA             Probability of False Alarm
%   PD              Probability of Detection
%   M               Number of samples collected
%
% Outputs:
%
%   xi              Required input signal-to-noise ratio [dB]
%
% Nicholas O'Donoughue
% 1 July 2019

% Check for existence of Statistics & Machine Learning Toolbox
use_stat_toolbox = license('test','Statistics_Toolbox');
   % If TRUE, then built-in functions will be used.
   % If FALSE, then custom-builts replacements in the utils namespace will
   % be used.
   
% Find the threshold
if use_stat_toolbox
    eta = chi2inv(1-PFA,2*M);
else
    eta = utils.chi2inv(1-PFA,2*M);
end

xi = zeros(size(eta));
for ii=1:numel(eta)
    thisEta = eta(ii);
    if numel(M)>1
        thisM = M(ii);
    else
        thisM = M;
    end
    
    if numel(PD)>1
        thisPD = PD(ii);
    else
        thisPD = PD;
    end
    
    % Set up function for probability of detection
    if use_stat_toolbox
        pd_fun = @(x) 1-ncx2cdf(thisEta,2*thisM,2*thisM*10.^(x/10)); % Xi is in dB
    else
        pd_fun = @(x) 1-utils.ncx2cdf(thisEta,2*thisM,2*thisM*10.^(x/10));
    end
    err_fun = @(x) pd_fun(x) - thisPD;
    
    % Initial Search Value
    thisXi = 0;             % Start at 0 dB
    err = err_fun(thisXi);  % Compute the difference between PD and desired PD
    
    % Initialize Search Parameters
    err_tol = .0001;    % Desired PD error tolerance
    maxIter = 1000;     % Maximum number of iterations
    idx = 0;            % Current iteration number
    maxStep=.5;          % Maximum Step Size [dB] - to prevent badly scaled
    % results when PD is near 0 or 1
    
    % Perform optimization
    while abs(err) > err_tol && idx < maxIter
        % Compute derivative at the current test point
        dxi = .01; % dB
        y0 = err_fun(thisXi-dxi);
        y1 = err_fun(thisXi+dxi);
        df = (y1-y0)/(2*dxi);
        
        % Newton-Rhapson Step Size
        step = -err/df;
        
        % Ensure that the step size is not overly large
        if abs(step)>maxStep
            step = sign(step)*maxStep;
        end
        
        % Iterate the Newton Approximation
        thisXi = thisXi + step;
        err = err_fun(thisXi);
        
        % Increment the iteration counter
        idx = idx + 1;
    end
    if idx>= maxIter
        warning('Computation finished before suitable tolerance achieved.  Error = %0.2f',abs(err));
    end
    xi(ii) = thisXi;
end