function [w, snrLoss] = makeTaper(N,taperType)
% [w, snrLoss] = makeTaper(N,taperType)
%
% Generate an amplitude taper of length N, according
% to the desired taperType, and optional set of parameters
%
% Inputs:
%
%   N           Length of the taper
%   taperType   String describing the type of taper desired.
%               Supported options are: "uniform", "cosine", "hanning",
%               "hamming", "bartlett", and "blackman-harris"
%
% Outputs:
%
%   w           Set of amplitude weights [0-1]
%   snrLoss     SNR Loss of peak return, w.r.t. uniform taper
%
% For discussion of these, and many other windows, see the Wikipedia page:
%   https://en.wikipedia.org/wiki/Window_function/
%
% Nicholas O'Donoughue
% 1 July 2019

switch lower(taperType)
    case 'uniform'
        w = ones(1,N);
        
    case 'cosine'
        i = 0:N-1;
        n = i - (N-1)/2;
        w = sin(pi/(2*N))*cos(pi*n/N);
        
    case 'hann'
        i = 0:N-1;
        n = i-(N-1)/2;
        w = cos(pi*n/N).^2;
        
    case 'hamming'
        i = 0:N-1;
        n = i - (N-1)/2;
        a = .54;%25/46;
        w = a + (1-a)*cos(2*pi*n/N);
        
    case 'blackman-harris'
        i = 0:N-1;
        n = i - (N-1)/2;
        w = 0.42 ...
            + 0.5*cos(2*pi*n/N) ...
            + 0.08*cos(4*pi*n/N);

    otherwise
        error('Invalid input; taper type %s not recognized.',taperType);
        
end

% Set peak to 1
w = w/max(abs(w(:)));

% Compute SNR Loss, rounded to the nearest hundredth of a dB
snrLoss = round(10*log10(sum(abs(w)/N)),2);