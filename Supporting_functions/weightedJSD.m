<<<<<<< HEAD
function jsd = weightedJSD(A,B,W)
% weightedJSD   Weighted Jensen–Shannon divergence that ignores NaNs
%
%   jsd = weightedJSD(A,B)          % uniform weights, ignore NaNs
%   jsd = weightedJSD(A,B,W)        % per-pixel weights (same size)
%
% OUTPUT
%   jsd   – scalar in [0 , ln(2)] ; NaN if < 2 valid pixels
%
% NOTE
%   Only the finite, non-NaN, positive-weight pixels contribute.
%   A and B are first normalized to probabilities over the valid set.

% ---------- defaults ----------------------------------------------------
if nargin<3 || isempty(W),  W = ones(size(A));  end
assert(isequal(size(A),size(B),size(W)), ...
       'A, B, and W must be the same size.');

A = double(A);  B = double(B);  W = double(W);

% ---------- validity mask ----------------------------------------------
mask =  isfinite(A) & isfinite(B) & isfinite(W) & (W>0);

if nnz(mask) < 2
    jsd = NaN;         % not enough data
    return
end

% ---------- apply mask & normalise -------------------------------------
A = A(mask);   B = B(mask);   W = W(mask);

W = W / sum(W);                       % normalization
A = A / sum(A);                       
B = B / sum(B);
M = 0.5*(A + B);

eps0 = 1e-12;                         % numerical guard 
term =  A .* log((A+eps0)./(M+eps0)) + ...
        B .* log((B+eps0)./(M+eps0));

jsd  = sum(W .* term);                
end
=======
function jsd = weightedJSD(A,B,W)
% weightedJSD   Weighted Jensen–Shannon divergence that ignores NaNs
%
%   jsd = weightedJSD(A,B)          % uniform weights, ignore NaNs
%   jsd = weightedJSD(A,B,W)        % per-pixel weights (same size)
%
% OUTPUT
%   jsd   – scalar in [0 , ln(2)] ; NaN if < 2 valid pixels
%
% NOTE
%   Only the finite, non-NaN, positive-weight pixels contribute.
%   A and B are first normalized to probabilities over the valid set.

% ---------- defaults ----------------------------------------------------
if nargin<3 || isempty(W),  W = ones(size(A));  end
assert(isequal(size(A),size(B),size(W)), ...
       'A, B, and W must be the same size.');

A = double(A);  B = double(B);  W = double(W);

% ---------- validity mask ----------------------------------------------
mask =  isfinite(A) & isfinite(B) & isfinite(W) & (W>0);

if nnz(mask) < 2
    jsd = NaN;         % not enough data
    return
end

% ---------- apply mask & normalise -------------------------------------
A = A(mask);   B = B(mask);   W = W(mask);

W = W / sum(W);                       % normalization
A = A / sum(A);                       
B = B / sum(B);
M = 0.5*(A + B);

eps0 = 1e-12;                         % numerical guard 
term =  A .* log((A+eps0)./(M+eps0)) + ...
        B .* log((B+eps0)./(M+eps0));

jsd  = sum(W .* term);                
end
>>>>>>> 5e4cf7b (Initial commit)
