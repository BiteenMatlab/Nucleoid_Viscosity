function [pSimu,pFit,simuResult,expResult,fSimu,fExp] = In_Nuc_fit( ...
            StepInNucSimu, StepInNuc, popNucRatio, ...
            slack_d, slack_acc, xscale, fig_plot)
%--------------------------------------------------------------------------
% In_Nuc_fit
% Fits step-size distributions obtained inside nucleoids with one- and two-
% population models.  First the simulated data are fit (1-pop) to obtain
% the parameters that later seed a 2-population fit of the experimental
% data.
%
% INPUTS
%   StepInNucSimu : N×3 numeric      – simulated trajectories; step sizes
%                                      must be in column 3
%   StepInNuc     : M×3 numeric      – experimental trajectories;
%   popNucRatio   : scalar [0-1]     – expected fraction of steps from nucleoid 
%
% Optional inputs :
%   slack_d       : double, default 0.15  – fractional bounds for Diffusion coefficient
%   slack_acc     : double, default 0.075 – fractional bounds for pop nuc ratio
%   xscale        : vector, default 0:20:1000 – histogram bin edges
%   fig_plot      : ["on","on"] | ["off","off"] | … – plot switches
%
% OUTPUTS
%   pSimu         : 1×2 fitted parameters for the simulated data
%   pFit          : 1×5 fitted parameters for the experimental data
%   simuResult    : struct from lsqcurvefit (resnorm, exitflag, output)
%   expResult     : experimental fit
%   fSimu/fExp    : figure handles (empty if plotting disabled)
%--------------------------------------------------------------------------
%
% Author: Xiaofeng Dai
% Date: 06/30/2025

%% ------------------------ 1. Handle optional inputs  --------------------
% Provide defaults then overwrite if the user supplied a value.
if nargin < 4 || isempty(slack_d),   slack_d   = 0.15;            end
if nargin < 5 || isempty(slack_acc), slack_acc = 0.075;           end
if nargin < 6 || isempty(xscale),    xscale    = 0:20:1000;       end
if nargin < 7 || isempty(fig_plot),  fig_plot  = ["on" "on"];     end

%% ------------------------ 2. Defensive checks ---------------------------
% Basic dimensional checks
assert(size(StepInNucSimu,2) >= 3, 'StepInNucSimu must have ≥3 columns');
assert(size(StepInNuc,    2) >= 3, 'StepInNuc     must have ≥3 columns');

% popNucRatio in [0,1]
assert(popNucRatio >= 0 && popNucRatio <= 1, ...
      'popNucRatio must be between 0 and 1');

% x-axis bins must be ascending
assert(all(diff(xscale) > 0), 'xscale has to be strictly increasing');

% Make fig_plot length-2 string array
if ischar(fig_plot) || isstring(fig_plot)
    fig_plot = repmat(string(fig_plot),1,2);
end

%% ------------------------ 3. PREP: common variables ---------------------
% Vectorised bin-centre calculation (avoids the loop)
binCenters = (xscale(1:end-1) + xscale(2:end))/2;

%% ======================= 4. ONE-POP fit (simulation) ====================
% Histogram (pdf normalised)
simuPDF = histcounts(StepInNucSimu(:,3), xscale, 'Normalization','pdf');

% Assemble data matrix expected by lsqcurvefit
xSimu = binCenters(:);
ySimu = simuPDF(:);

% 1-pop functional form
simuFun = @(p,x) (2*x/p(1)).*exp(-(x.^2)/p(1)) + p(2)*x;

% Initial guess & bounds
p0_simu = [4.8e4, 0];
lb_simu = [1e4,  0];
ub_simu = [1e5,  0];    % 2nd parameter is forced to 0 here

% Solver options (tight tolerances)
opt = optimoptions('lsqcurvefit', 'Algorithm','trust-region-reflective', ...
                   'Display','off', 'MaxIterations',400, ...
                   'OptimalityTolerance',1e-9, 'FunctionTolerance',1e-9, ...
                   'StepTolerance',1e-13);

[pSimu, resnorm,~,exitflag,out] = lsqcurvefit(simuFun, p0_simu, ...
                                             xSimu, ySimu, lb_simu, ub_simu, opt);

% Package diagnostic information
simuResult = struct('resnorm',resnorm,'exitflag',exitflag,'output',out);

% Plot simulation fit if requested
if strcmpi(fig_plot(1),'on')
    fSimu = figure('Name','Simulated displacement fitting');
    histogram(StepInNucSimu(:,3), xscale, 'Normalization','pdf');
    hold on
    plot(xSimu, simuFun(pSimu,xSimu), 'r--','LineWidth',1.5)
    xlabel('Displacement'); ylabel('PDF'); legend('Histogram','Fit');
    title({'Displacement from Simulation';'1-population fit'})
else
    fSimu = [];
end


%% ======================= 5. TWO-POP fit (experiment) ====================
% Experimental histogram
expPDF = histcounts(StepInNuc(:,3), xscale, 'Normalization','pdf');
xExp   = binCenters(:);
yExp   = expPDF(:);

% 2-population model
expFun = @(p,x) ...
         ( (2*x/p(2)).*exp(-(x.^2)/p(2)) + p(3)*x )          * p(1) + ...
         ( (2*x/p(4)).*exp(-(x.^2)/p(4)) + p(5)*x ) .* (1-p(1));

% Initial guess and bounds
p0_exp = [popNucRatio, 1.4e4, 0,  pSimu(1), pSimu(2)];
lb_exp = [max(0, popNucRatio-slack_acc),  3e3, 0, ...
          pSimu(1)*(1-slack_d),           pSimu(2)*(1-slack_d)];
ub_exp = [min(1, popNucRatio+slack_acc), 3.2e4, 0, ...
          pSimu(1)*(1+slack_d),           pSimu(2)*(1+slack_d)];

[pFit, resnorm,~,exitflag,out] = lsqcurvefit(expFun, p0_exp, ...
                                     xExp, yExp, lb_exp, ub_exp, opt);

expResult = struct('resnorm',resnorm,'exitflag',exitflag,'output',out);

% Plot experimental fit if requested
if strcmpi(fig_plot(2),'on')
    fExp = figure('Name','Experimental displacement fitting');
    histogram(StepInNuc(:,3), xscale, 'Normalization','pdf','DisplayName','Histogram');
    hold on

    yPop1 = ( (2*xExp/pFit(2)).*exp(-(xExp.^2)/pFit(2)) + pFit(3)*xExp ) * pFit(1);
    yPop2 = ( (2*xExp/pFit(4)).*exp(-(xExp.^2)/pFit(4)) + pFit(5)*xExp ) * (1-pFit(1));
    yRecon = yPop1 + yPop2;
    yErr   = yExp  - yRecon;

    plot(xExp, yPop1,'y--','LineWidth',1.5,'DisplayName','Pop 1');
    plot(xExp, yPop2,'r--','LineWidth',1.5,'DisplayName','Pop 2');
    plot(xExp, yRecon,'m-','LineWidth',1,'DisplayName','Pop 1+2');
    plot(xExp, yErr,'b--','LineWidth',1.5,'DisplayName','Residual');
    title({'Displacement from Experiment';'2-population fit'});
    xlabel('Step length'); ylabel('PDF');
    legend show
    ylim([1.5*min(yErr) 1.1*max(yRecon)])
else
    fExp = [];
end
end