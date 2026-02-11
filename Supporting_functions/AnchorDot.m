function [anchor_coord, slope] = AnchorDot(outline)
% AnchorDot Reduces an outline by keeping only anchor points along collinear segments.
%
% [anchor_coord, slope] = AnchorDot(outline) takes an N-by-2 array of [x,y]
% coordinates representing a boundary (or outline) and removes redundant points
% that lie on the same (approximately) straight line. Only the starting and ending
% points of collinear segments are retained. The function also returns the slopes
% between consecutive anchor points.
%
% Input:
% outline - An N-by-2 numeric array of [x,y] coordinates.
%
% Outputs:
% anchor_coord - An M-by-2 numeric array containing the anchor (reduced) coordinates.
% slope - An (M-1)-by-1 numeric vector containing slopes between consecutive anchors.
%
% Example:
% [anchor, segSlope] = AnchorDot(outline);
%
% Remarks:
% A tolerance is used when comparing slopes to decide if three consecutive points
% are collinear.
%
% Author: Xiaofeng Dai
% Date: 04/15/2025

%% 1. Input Validation
if ~isnumeric(outline) || size(outline,2) ~= 2
    error('Input "outline" must be an N-by-2 numeric array.');
end
N = size(outline,1);
if N < 2
    anchor_coord = outline;
    slope = [];
    return;
end

%% 2. Set Tolerance and Initialize Anchor Indices
tol = 1e-3;  % Tolerance for detecting differences in slope
anchorIdx = 1;  % Always include the first point

%% 3. Loop Through Points to Detect Change in Slope
% For each potential anchor point (from second to next-to-last),
% compute the slopes of segments [i-1,i] and [i,i+1] and compare.
for i = 2:(N-1)
    % Slope from point i-1 to i:
    delta1 = outline(i,:) - outline(i-1,:);
    if abs(delta1(1)) < 1e-10
        slope1 = sign(delta1(2)) * inf;
    else
        slope1 = delta1(2) / delta1(1);
    end

    % Slope from point i to i+1:
    delta2 = outline(i+1,:) - outline(i,:);
    if abs(delta2(1)) < 1e-10
        slope2 = sign(delta2(2)) * inf;
    else
        slope2 = delta2(2) / delta2(1);
    end

    % If the difference between slopes exceeds tolerance, point i
    % is not collinear with its neighbors.
    if abs(slope1 - slope2) > tol
        anchorIdx(end+1) = i;
    end
end

% Always include the last point.
anchorIdx(end+1) = N;

%% 4. Extract Anchor Coordinates and Compute Segment Slopes
anchor_coord = outline(anchorIdx, :);
M = length(anchorIdx);
slope = zeros(M-1, 1);
for j = 1:(M-1)
    diffVec = anchor_coord(j+1,:) - anchor_coord(j,:);
    if abs(diffVec(1)) < 1e-10
        slope(j) = sign(diffVec(2)) * inf;
    else
        slope(j) = diffVec(2) / diffVec(1);
    end
end
end
