<<<<<<< HEAD
function ciNuc = In_Nuc_fit_ci(stepInNuc, pFit, xscale, alpha)
%--------------------------------------------------------------------------
% In_Nuc_fit_ci
% Confidence interval for the “a_nuc” parameter (p_fit(2)) obtained in the
% two-population fit performed by In_Nuc_fit.m.
%
% INPUTS
%   stepInNuc : N×3 numeric   – experimental steps (col-3 is step length)
%   pFit      : 1×5 numeric   – parameter vector returned by In_Nuc_fit
%
% Optional
%   xscale    : 1×K+1 numeric – histogram bin edges    (default 0:20:1000)
%   alpha     : scalar        – significance level     (default 0.05 ⇒ 95 % CI)
%
% OUTPUT
%   ciNuc     : 2×1 numeric   – [lowerBound; upperBound] for a_nuc
%--------------------------------------------------------------------------
%
% Author: Xiaofeng Dai
% Date: 06/30/2025

%% ------------------------ 1. Defaults & validation ----------------------
if nargin < 3 || isempty(xscale), xscale = 0:20:1000; end
if nargin < 4 || isempty(alpha),  alpha  = 0.05;      end

% Defensive checks
assert(size(stepInNuc,2) >= 3, '`stepInNuc` must have ≥3 columns.');
assert(numel(pFit)       == 5, '`pFit` must have 5 elements.');
assert(all(diff(xscale)  > 0), '`xscale` must be strictly increasing.');
assert(alpha > 0 && alpha < 1, '`alpha` must be in the interval (0,1).');

%% ------------------------ 2. Build histogram ---------------------------
binCtr = (xscale(1:end-1) + xscale(2:end))/2;                         % vectorised
pdfVals = histcounts(stepInNuc(:,3), xscale,'Normalization','pdf')';

%% ------------------------ 3. Prepare model -----------------------------
weight  = pFit(1);           % fixed population weight
aOut    = pFit(4);           % fixed “outside” parameter

% The model has a single free parameter: a_nuc
modelFun = @(a,x) ( 2*x./a   .*exp(-(x.^2)./a)   ) * weight + ...
                  ( 2*x./aOut.*exp(-(x.^2)./aOut)) * (1-weight);

% Initial guess: best estimate from previous fit
a0 = pFit(2);

%% ------------------------ 4. Non-linear fit + CI ------------------------
% Tight-ish options to match the main fit routine
opts = statset('nlinfit');
opts.RobustWgtFun = [];      % same as default – can be changed if needed
opts.MaxIter      = 400;
opts.TolX         = 1e-12;
opts.TolFun       = 1e-12;

[beta, residual, ~, covB] = nlinfit(binCtr', pdfVals, modelFun, a0, opts);

ciNuc = nlparci(beta, residual, 'Covar', covB, 'Alpha', alpha);
ciNuc = ciNuc(:);            % return as 2×1 column vector
=======
function ciNuc = In_Nuc_fit_ci(stepInNuc, pFit, xscale, alpha)
%--------------------------------------------------------------------------
% In_Nuc_fit_ci
% Confidence interval for the “a_nuc” parameter (p_fit(2)) obtained in the
% two-population fit performed by In_Nuc_fit.m.
%
% INPUTS
%   stepInNuc : N×3 numeric   – experimental steps (col-3 is step length)
%   pFit      : 1×5 numeric   – parameter vector returned by In_Nuc_fit
%
% Optional
%   xscale    : 1×K+1 numeric – histogram bin edges    (default 0:20:1000)
%   alpha     : scalar        – significance level     (default 0.05 ⇒ 95 % CI)
%
% OUTPUT
%   ciNuc     : 2×1 numeric   – [lowerBound; upperBound] for a_nuc
%--------------------------------------------------------------------------
%
% Author: Xiaofeng Dai
% Date: 06/30/2025

%% ------------------------ 1. Defaults & validation ----------------------
if nargin < 3 || isempty(xscale), xscale = 0:20:1000; end
if nargin < 4 || isempty(alpha),  alpha  = 0.05;      end

% Defensive checks
assert(size(stepInNuc,2) >= 3, '`stepInNuc` must have ≥3 columns.');
assert(numel(pFit)       == 5, '`pFit` must have 5 elements.');
assert(all(diff(xscale)  > 0), '`xscale` must be strictly increasing.');
assert(alpha > 0 && alpha < 1, '`alpha` must be in the interval (0,1).');

%% ------------------------ 2. Build histogram ---------------------------
binCtr = (xscale(1:end-1) + xscale(2:end))/2;                         % vectorised
pdfVals = histcounts(stepInNuc(:,3), xscale,'Normalization','pdf')';

%% ------------------------ 3. Prepare model -----------------------------
weight  = pFit(1);           % fixed population weight
aOut    = pFit(4);           % fixed “outside” parameter

% The model has a single free parameter: a_nuc
modelFun = @(a,x) ( 2*x./a   .*exp(-(x.^2)./a)   ) * weight + ...
                  ( 2*x./aOut.*exp(-(x.^2)./aOut)) * (1-weight);

% Initial guess: best estimate from previous fit
a0 = pFit(2);

%% ------------------------ 4. Non-linear fit + CI ------------------------
% Tight-ish options to match the main fit routine
opts = statset('nlinfit');
opts.RobustWgtFun = [];      % same as default – can be changed if needed
opts.MaxIter      = 400;
opts.TolX         = 1e-12;
opts.TolFun       = 1e-12;

[beta, residual, ~, covB] = nlinfit(binCtr', pdfVals, modelFun, a0, opts);

ciNuc = nlparci(beta, residual, 'Covar', covB, 'Alpha', alpha);
ciNuc = ciNuc(:);            % return as 2×1 column vector
>>>>>>> 5e4cf7b (Initial commit)
end