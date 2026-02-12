<<<<<<< HEAD
function [score,Smap] = weightedSSIM(A,B,W,radius)
% weightedSSIM  –  NaN-safe weighted SSIM for two 2-D matrices
%
%   score            = weightedSSIM(A,B,W)
%   [score,Smap]     = weightedSSIM(A,B,W,radius)
%
% INPUTS
%   A, B   : M×N matrices.   May contain NaNs.
%   W      : M×N non-negative weights; [] → uniform.
%   radius : Gaussian window radius (default 2.5 → 6×6 window).
%
% OUTPUTS
%   score  : global weighted SSIM  (-1 ~ 1; 1 = identical).
%   Smap   : local SSIM map   (optional).

% ---------- defaults ----------------------------------------------------
if nargin<3 || isempty(W),  W = ones(size(A),'like',A);  end
if nargin<4 || isempty(radius), radius = 2.5;              end
assert(isequal(size(A),size(B),size(W)), ...
      'A, B, and W must have identical sizes.');

A = double(A);  B = double(B);  W = double(W);

% ---------- build validity mask & zero-out bad pixels -------------------
mask =  isfinite(A) & isfinite(B) & isfinite(W) & (W>0);
if nnz(mask) < 2
    score = NaN;              % not enough valid information
    if nargout>1, Smap = NaN(size(A)); end
    return
end

A(~mask) = 0;   B(~mask) = 0;   W(~mask) = 0;           % neutralise NaNs
W = W / sum(W(:));                                      % ΣW = 1

% ---------- dynamic range using only valid pixels -----------------------
L = max([A(mask); B(mask)]) - min([A(mask); B(mask)]);
if L==0, L = max(std([A(mask); B(mask)]), eps); end

% ---------- local SSIM map ---------------------------------------------
[~,Smap] = ssim(A,B, ...
                'DynamicRange',L, ...
                'Radius',      radius, ...
                'Exponents',   [1 1 1]);

% ---------- smooth & re-normalise weight mask --------------------------
gwin = fspecial('gaussian', 2*radius+1, 1.5);   % same σ as SSIM default
Wsm  = conv2(W, gwin, 'same');
Wsm(~mask) = 0;                                 % keep bad pixels at zero
Wsm  = Wsm / max(sum(Wsm(:)), eps);             % ΣWsm = 1

% ---------- global weighted SSIM ---------------------------------------
score = sum(Wsm .* Smap, "all");
end
=======
function [score,Smap] = weightedSSIM(A,B,W,radius)
% weightedSSIM  –  NaN-safe weighted SSIM for two 2-D matrices
%
%   score            = weightedSSIM(A,B,W)
%   [score,Smap]     = weightedSSIM(A,B,W,radius)
%
% INPUTS
%   A, B   : M×N matrices.   May contain NaNs.
%   W      : M×N non-negative weights; [] → uniform.
%   radius : Gaussian window radius (default 2.5 → 6×6 window).
%
% OUTPUTS
%   score  : global weighted SSIM  (-1 ~ 1; 1 = identical).
%   Smap   : local SSIM map   (optional).

% ---------- defaults ----------------------------------------------------
if nargin<3 || isempty(W),  W = ones(size(A),'like',A);  end
if nargin<4 || isempty(radius), radius = 2.5;              end
assert(isequal(size(A),size(B),size(W)), ...
      'A, B, and W must have identical sizes.');

A = double(A);  B = double(B);  W = double(W);

% ---------- build validity mask & zero-out bad pixels -------------------
mask =  isfinite(A) & isfinite(B) & isfinite(W) & (W>0);
if nnz(mask) < 2
    score = NaN;              % not enough valid information
    if nargout>1, Smap = NaN(size(A)); end
    return
end

A(~mask) = 0;   B(~mask) = 0;   W(~mask) = 0;           % neutralise NaNs
W = W / sum(W(:));                                      % ΣW = 1

% ---------- dynamic range using only valid pixels -----------------------
L = max([A(mask); B(mask)]) - min([A(mask); B(mask)]);
if L==0, L = max(std([A(mask); B(mask)]), eps); end

% ---------- local SSIM map ---------------------------------------------
[~,Smap] = ssim(A,B, ...
                'DynamicRange',L, ...
                'Radius',      radius, ...
                'Exponents',   [1 1 1]);

% ---------- smooth & re-normalise weight mask --------------------------
gwin = fspecial('gaussian', 2*radius+1, 1.5);   % same σ as SSIM default
Wsm  = conv2(W, gwin, 'same');
Wsm(~mask) = 0;                                 % keep bad pixels at zero
Wsm  = Wsm / max(sum(Wsm(:)), eps);             % ΣWsm = 1

% ---------- global weighted SSIM ---------------------------------------
score = sum(Wsm .* Smap, "all");
end
>>>>>>> 5e4cf7b (Initial commit)
