function [pCytoFit, pCytoCI, cytoResult, hFig] = In_Cyto_fit( ...
                  rawSpatialDisp, spMap, coorInCyto, ...
                  asptRatio, numPixY, xscale, fig_plot, plot_err)
%--------------------------------------------------------------------------
% In_Cyto_fit
% One-population fit of step-size distribution inside the cytoplasm.
%
% INPUTS
%   rawSpatialDisp : N×3 numeric   – [x, y, stepLength] for every step
%   spMap          : spatial Displacement Map  
%   coorInCyto     : K×2  int      – list of [x y] pixel indices inside cyto
%   asptRatio      : scalar        – aspect ratio of cell (x-scale factor)
%
% Optional
%   numPixY        : scalar        – #pixels along y  (default 20)
%   xscale         : vector        – histogram bin-edges (default 0:20:1000)
%   fig_plot       : 'on'/'off'    – show main plot     (default 'on')
%   plot_err       : 'on'/'off'    – overlay residual   (default 'off')
%
% OUTPUTS
%   pCytoFit   : 1×2 double   – fitted parameters
%   pCytoCI    : 2×1 double   – 95 % CI for a_cyto   (change alpha inside)
%   cytoResult : struct       – diagnostics from lsqcurvefit
%   hFig       : figure handle (empty if fig_plot = 'off')
%--------------------------------------------------------------------------
%
% Author: Xiaofeng Dai
% Date: 06/30/2025

%% ------------------------ 1. Defaults & validation ----------------------
if nargin < 5 || isempty(numPixY), numPixY = 20;          end
if nargin < 6 || isempty(xscale),  xscale  = 0:20:1000;   end
if nargin < 7 || isempty(fig_plot),fig_plot = 'on';       end
if nargin < 8 || isempty(plot_err),plot_err = 'off';      end

% Basic dimension checks
assert(size(rawSpatialDisp,2) == 3, ...
       '`rawSpatialDisp` must be N×3: [x y stepLength].');
assert(ismatrix(spMap) && ~isempty(spMap), ...
       '`spMap` must be a 2-D numeric matrix.');
assert(size(coorInCyto,2) == 2, '`coorInCyto` must be K×2.');
assert(all(diff(xscale) > 0),   '`xscale` must be strictly increasing.');

%% ------------------------ 2. Effective-pixel coordinates ----------------
% Non-NaN locations in spMap are “valid” pixels
[row, col] = find(~isnan(spMap));
coorEff    = [col, row];               % [x y] pairs (x = col, y = row)

% Cytoplasm coordinates that overlap with valid pixels
cytoCoord = coorInCyto(ismember(coorInCyto, coorEff, 'rows'), :);

%% ------------------------ 3. Extract steps that fall inside cytoplasm ---
stepInCyto = Extract_Step_InRegion(rawSpatialDisp, cytoCoord, ...
                                   numPixY, asptRatio);

%% ------------------------ 4. Build histogram (PDF-normalised) ----------
binCtr  = (xscale(1:end-1) + xscale(2:end)) / 2;          % vectorised
cytoPDF = histcounts(stepInCyto(:,3), xscale, 'Normalization','pdf')';

xFit = binCtr(:);
yFit = cytoPDF(:);

%% ------------------------ 5. One-population least-squares fit ----------
fitFun = @(p,x) (2*x./p(1)).*exp(-(x.^2)./p(1)) + p(2)*x;  
p0    = [2.5e4, 0];
lb    = [5e3 , 0];
ub    = [1e5 , 5e-8];

opt = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective', ...
                   'Display','off','MaxIter',400, ...
                   'OptimalityTolerance',1e-10, ...
                   'FunctionTolerance',1e-10, 'StepTolerance',1e-14);

[pCytoFit, resnorm, ~, exitflag, out] = ...
            lsqcurvefit(fitFun, p0, xFit, yFit, lb, ub, opt);

cytoResult = struct('resnorm',resnorm,'exitflag',exitflag,'output',out);

%% ------------------------ 6. Confidence interval for a_cyto -------------
singlePop = @(a,x) (2*x./a).*exp(-(x.^2)./a);    % slope term fixed at 0
alphaCI   = 0.05;                                % 95 % CI
[beta, resid, ~, covB] = nlinfit(xFit, yFit, singlePop, pCytoFit(1));
pCytoCI   = nlparci(beta, resid, 'Covar', covB, 'Alpha', alphaCI);
pCytoCI   = pCytoCI(:);                          % column vector [low; high]

%% ------------------------ 7. Plot (if requested) ------------------------
if strcmpi(fig_plot,'on')
    hFig = figure('Name','Cytoplasm – 1-population fit');
    histogram(stepInCyto(:,3), xscale, 'Normalization','pdf'); hold on
    plot(xFit, fitFun(pCytoFit,xFit), 'r--','LineWidth',1.5);

    if strcmpi(plot_err,'on')
        plot(xFit, yFit - fitFun(pCytoFit,xFit), ...
             'b--','LineWidth',1.2, 'DisplayName','Residual');
        ylim([1.1*min(yFit - fitFun(pCytoFit,xFit)), ...
              1.1*max(max(fitFun(pCytoFit,xFit)), max(yFit))]);
    end
    xlabel('Step Size'); ylabel('PDF');
    legend('Data','Fit','Location','best');
    title('Displacement from Cytoplasm')
else
    hFig = [];
end
end