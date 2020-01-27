function [gammaOx,gammaW] = gasLossCoeff(f,p,e,T)
% [gammaOx,gammaW] = gasLossCoeff(f,P,e,T)
%
% Inputs:
%   f       - Propagation Frequencies [Hz]
%   P       - Dry Air Pressure [hPa]
%   e       - Water Vapor Partial Pressure [hPa]
%   T       - Temperature [K]
%
% Outputs:
%   gammaOx - Gas loss coefficient due to oxygen [dB/km]
%   gammaW  - Gas loss coefficient due to water vapor [dB/km]
%
% Implement the atmospheric loss coefficients from Annex 1 of ITU-R
% P.676-11 (12/2017)
%
% If array inputs are specified, then array results are given for
% alphaO and alphaW.
%
% 1 July 2019
% Nicholas O'Donoughue

if all(f(:)==0 | isnan(f(:)))
    gammaOx=0;
    gammaW=0;
    return;
end

% Account for vector inputs
if numel(f) > 1 || numel(p) > 1 || numel(e) > 1 || numel(T) > 1
    % At least one input is non-scalar, find the dimensions
    dimsF = size(f);
    dimsP = size(p);
    dimsE = size(e);
    dimsT = size(T);
    
    % Check for edge case, all inputs either scalar or equal dim
    maxDim = max([dimsF;dimsP;dimsE;dimsT],[],1); % Find max dimension
    if (all(dimsF==maxDim) || numel(f)==1) && ...
       (all(dimsP==maxDim) || numel(p)==1) && ...
       (all(dimsE==maxDim) || numel(e)==1) && ...
       (all(dimsT==maxDim) || numel(T)==1)
   
        % Reshape all inputs to a row vector, since spectroscopic data
        % will be a column vector
        f = f(:).';
        p = p(:).';
        e = e(:).';
        T = T(:).';
    else
        % Some inputs are differently sized, expand all non-zero along
        % dims
        if numel(f) > 1
            % Check that all non-matching dims are scalar
            nonMatchDim = dimsF~=maxDim;
            assert(all(dimsF(nonMatchDim)==1),'Error: All non-singleton dimensions of inputs must match.');
            
            % Find expansion
            repDim = ones(size(dimsF));
            repDim(nonMatchDim) = maxDim(nonMatchDim);
            f = repmat(f,repDim);
            f = f(:).';
        end
        
        if numel(p) > 1
            % Check that all non-matching dims are scalar
            nonMatchDim = dimsP~=maxDim;
            assert(all(dimsP(nonMatchDim)==1),'Error: All non-singleton dimensions of inputs must match.');
            
            % Find expansion
            repDim = ones(size(dimsP));
            repDim(nonMatchDim) = maxDim(nonMatchDim);
            p = repmat(p,repDim);
            p = p(:).';
        end
        
        if numel(e) > 1
            % Check that all non-matching dims are scalar
            nonMatchDim = dimsE~=maxDim;
            assert(all(dimsE(nonMatchDim)==1),'Error: All non-singleton dimensions of inputs must match.');
            
            % Find expansion
            repDim = ones(size(dimsE));
            repDim(nonMatchDim) = maxDim(nonMatchDim);
            e = repmat(e,repDim);
            e = e(:).';
        end
        
        if numel(T) > 1
            % Check that all non-matching dims are scalar
            nonMatchDim = dimsT~=maxDim;
            assert(all(dimsT(nonMatchDim)==1),'Error: All non-singleton dimensions of inputs must match.');
            
            % Find expansion
            repDim = ones(size(dimsT));
            repDim(nonMatchDim) = maxDim(nonMatchDim);
            T = repmat(T,repDim);
            T = T(:).';
        end
    end
    outputDims = maxDim;
else
    outputDims = [1 1];
end
% At this point, all of the inputs are either scalar, or 1xM

% Read in the spectroscopic tables (Tables 1 and 2 of Annex 1)
% All table data will be Nx1 for the N spectroscopic lines of
% each table.
[fox,a1,a2,a3,a4,a5,a6] = atm.makeSpectroscopicTableOxygen();
% [fox,a1,a2,a3,a4,a5,a6]=textread('+atm/spectroscopicTableOxygen.dat','%f,%f,%f,%f,%f,%f,%f','headerlines',1);
[fw,b1,b2,b3,b4,b5,b6] = atm.makeSpectroscopicTableWater();
% [fw,b1,b2,b3,b4,b5,b6]=textread('+atm/spectroscopicTableWater.dat','%f,%f,%f,%f,%f,%f,%f','headerlines',1);


%% Compute the dry continuum due to pressure-induced Nitrogen absorption and the Debye spectrum (eq 8)
f0 = f/1e9; % Convert freq from Hz to GHz
th = 300./T;
d = 5.6e-4*(p+e).*th.^.8;
ND = f0.*p.*th.^2.*(6.14e-5./(d.*(1+(f0./d).^2))+(1.4e-12*p.*th.^1.5)./(1+1.9e-5*f0.^1.5));

%% Compute the strength of the i-th water/o2 vapor line (eq 3)
Sox = (a1.*1e-7.*p.*th.^3).*exp(a2.*(1-th));
Sw = (b1.*1e-1.*e.*th.^3.5).*exp(b2.*(1-th));

%% Compute the line shape factor for each
%  Correction factor due to interference effects in oxygen lines (eq 7)
dox = (a5+a6.*th).*1e-4.*(p+e).*th.^.8;
dw = 0;

% spectroscopic line width (eq 6a)
dfox = (a3*1e-4).*(p.*th.^(.8-a4)+1.1*e.*th);
dfw = (b3*1e-4).*(p.*th.^(b4)+b5.*e.*th.^b6);

% modify spectroscopic line width to account for Zeeman splitting of oxygen
% lines and Doppler broadening of water vapour lines (eq 6b)
dfox_sq = dfox.^2+2.25e-6;
dfox = sqrt(dfox_sq);
dfw = .535.*dfw+sqrt(.217*dfw.^2+(2.1316e-12.*fw.^2)./th);

% Compute line shape factor
delta_fox = fox-f0;
sum_fox = fox+f0;
delta_fw = fw-f0;
sum_fw = fw+f0;
Fox = f0./fox.*( ((dfox-dox.*delta_fox)./(delta_fox.^2+dfox.^2)) + ((dfox-dox.*sum_fox)./(sum_fox.^2+dfox_sq)));
Fw = f0./fw.*( ((dfw-dw.*delta_fw)./(delta_fw.^2+dfw.^2)) + ((dfw-dw.*sum_fw)./(sum_fw.^2+dfw.^2)));

%% Compute complex refractivities
Nox = sum(Sox.*Fox,1)+ND;
Nw = sum(Sw.*Fw,1);

gammaOx = reshape(.1820*f0.*Nox,outputDims);
gammaW = reshape(.1820*f0.*Nw,outputDims);

%% Handle all freqs < 1 GHz
lowFreq = f0 < 1;
gammaOx(lowFreq) = 0;
gammaW(lowFreq) = 0;
