<<<<<<< HEAD
function [stepInRegion, stepOut] = Extract_Step_InRegion( ...
                        rawSpatialDisp, coord, numPixY, asptRatio)
%--------------------------------------------------------------------------
% Extract_Step_InRegion
% Returns the subset of steps whose (x,y) pixel indices fall in `coord`.
%
% OUTPUTS
%   stepInRegion : M×3 [xIdx yIdx stepLength]
%   stepOut      : R×3 steps outside the region (not used by caller)
%--------------------------------------------------------------------------

% Pixel size (fraction of cell length) along each axis
xPix = 1 / (asptRatio * numPixY);
yPix = 1 /  numPixY;

% Vectorised mapping from real-valued coords → pixel indices
xIdx = floor((rawSpatialDisp(:,1) + 0.5) / xPix) + 1;
yIdx = floor((rawSpatialDisp(:,2) + 0.5) / yPix) + 1;

idxPair = [xIdx, yIdx];

inMask  = ismember(idxPair, coord, 'rows');   % logical mask
stepInRegion = [idxPair(inMask ,:) , rawSpatialDisp(inMask ,3)];
stepOut      = [idxPair(~inMask,:) , rawSpatialDisp(~inMask,3)];
=======
function [stepInRegion, stepOut] = Extract_Step_InRegion( ...
                        rawSpatialDisp, coord, numPixY, asptRatio)
%--------------------------------------------------------------------------
% Extract_Step_InRegion
% Returns the subset of steps whose (x,y) pixel indices fall in `coord`.
%
% OUTPUTS
%   stepInRegion : M×3 [xIdx yIdx stepLength]
%   stepOut      : R×3 steps outside the region (not used by caller)
%--------------------------------------------------------------------------

% Pixel size (fraction of cell length) along each axis
xPix = 1 / (asptRatio * numPixY);
yPix = 1 /  numPixY;

% Vectorised mapping from real-valued coords → pixel indices
xIdx = floor((rawSpatialDisp(:,1) + 0.5) / xPix) + 1;
yIdx = floor((rawSpatialDisp(:,2) + 0.5) / yPix) + 1;

idxPair = [xIdx, yIdx];

inMask  = ismember(idxPair, coord, 'rows');   % logical mask
stepInRegion = [idxPair(inMask ,:) , rawSpatialDisp(inMask ,3)];
stepOut      = [idxPair(~inMask,:) , rawSpatialDisp(~inMask,3)];
>>>>>>> 5e4cf7b (Initial commit)
end