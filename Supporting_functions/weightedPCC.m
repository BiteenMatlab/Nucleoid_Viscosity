function r = weightedPCC(A,B,W)
% weightedPCC   Weighted Pearson correlation that skips NaNs
%
%   r = weightedPCC(A,B)          % uniform weights, ignore NaNs
%   r = weightedPCC(A,B,W)        % per-pixel weights (same size as A,B)
%
% OUTPUT
%   r  – scalar in [-1, 1] ;
%


    % ---------- defaults & checks ----------
    if nargin < 3 || isempty(W),  W = ones(size(A));  end
    assert(isequal(size(A),size(B),size(W)), ...
          'A, B, W must be the same size.');

    % ---------- build validity mask ----------
    mask =  ~isnan(A) & ~isnan(B) & ~isfinite(A)==0 & ~isfinite(B)==0 ...
            & ~isnan(W) & (W > 0);

    if nnz(mask) < 2
        r = NaN;     % not enough valid points
        return
    end

    % ---------- apply mask & reshape ----------
    A = A(mask);   B = B(mask);   W = W(mask);

    % ---------- normalise weights ----------
    W = W / sum(W);

    % ---------- weighted means ----------
    muA = sum(W .* A);
    muB = sum(W .* B);

    % ---------- numerator & denominator ----------
    num = sum(W .* (A - muA) .* (B - muB));
    den = sqrt( sum(W .* (A - muA).^2) * sum(W .* (B - muB).^2) );

    % ---------- correlation ----------
    r = num / den;
end
