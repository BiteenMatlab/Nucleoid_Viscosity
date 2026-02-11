function [nucVolumeRatio, cmpAccHeatMap] = Nuc_acc(genDot, accRange, R, L, anchorModi, scalingWidth, scalingLength, numPixY, aspectRatio, numDot, figPlot)
% Nuc_acc simulates dot localization within a cell and its nucleoid,
% calculates the nucleoid volume ratio, and generates accessibility heat maps.
%
% Inputs:
%   genDot      - Total number of generated dots.
%   accRange    - Vector of accessibility range factors.
%   R           - Radius parameter for the cell.
%   L           - Length offset for the cell.
%   anchorModi  - Matrix defining nucleoid boundary anchor points (minimum 2 rows).
%   scalingWidth- Scaling factor for the x-dimension.
%   scalingLength- Scaling factor for the y-dimension.
%   numPixY     - Number of pixels along the y axis for the heatmap.
%   aspectRatio - Aspect ratio for the heatmap.
%   numDot      - Number of dots to use for the heat map.
%   figPlot     - (Optional) 'on' to display plots, 'off' to suppress; default 'off'.
%
% Outputs:
%   nucVolumeRatio  - Estimated nucleoid volume ratio.
%   cmpAccHeatMap   - Cell array of comparison accessibility heat maps.
%
% Example usage:
%   [volRatio, heatMaps] = Nuc_acc(1e6, (0:0.1:1), 5, 3, anchorPts, 6, 12.5, 20, 2, 4e5, 'on');
% Author: Xiaofeng Dai
% Date: 05/17/2025

% Check number of inputs and provide default for figPlot if necessary
if nargin < 10
    error('Nuc_acc requires at least 10 arguments.');
elseif nargin < 11
    figPlot = 'off';
end

% ----- Defensive Checks -----
if ~isscalar(genDot) || genDot <= 0
    error('genDot must be a positive scalar.');
end
if ~isvector(accRange) || isempty(accRange)
    error('accRange must be a non-empty vector.');
end
if ~isscalar(R) || R <= 0
    error('R must be a positive scalar.');
end
if ~isscalar(L) || L < 0
    error('L must be a nonnegative scalar.');
end
if size(anchorModi,1) < 2
    error('anchorModi must contain at least two rows to define a boundary.');
end
if ~isscalar(scalingWidth) || scalingWidth <= 0 || ~isscalar(scalingLength) || scalingLength <= 0
    error('Scaling factors must be positive scalars.');
end
if ~isscalar(numPixY) || numPixY <= 0
    error('numPixY must be a positive scalar.');
end
if ~isscalar(numDot) || numDot <= 0
    error('numDot must be a positive scalar.');
end
if genDot <= numDot
    error('The number of generated dots (genDot) must exceed the number of dots used for the heat map (numDot).');
end

% ----- Compute Nucleoid Volume Ratio and Dot Coordinates -----
[nucVolumeRatio, dotCoordinates] = Cal_nuc_vol(genDot, R, L, anchorModi);

% Optionally display the nucleoid volume ratio in the command window.
if strcmpi(figPlot, 'on')
    fprintf('Nucleoid Volume Ratio = %.4f (N = %d)\n', nucVolumeRatio, genDot);
end

% Separate dot coordinates by region
%   Region flag 0     : Cytoplasm 
%   Region flag 2     : Nucleoid
dotsCytoplasm = dotCoordinates(dotCoordinates(:,4) == 0, :);
dotsNucleoid  = dotCoordinates(dotCoordinates(:,4) == 2, :);

% Calculate the number of dots to randomly sample from each region
numDotsNuc   = floor(accRange .* nucVolumeRatio * numDot ./ (1 - nucVolumeRatio + accRange .* nucVolumeRatio));
numDotsCyto  = floor((1 - nucVolumeRatio) * numDot ./ (1 - nucVolumeRatio + accRange .* nucVolumeRatio));

% Preallocate cell array for the coordinate subsets at each accRange value.
coordForAccessibility = cell(length(accRange), 1);
for idx = 1:length(accRange)
    % Defensive check: if requested dots exceed available dots, use the available maximum.
    nucRequested = numDotsNuc(idx);
    cytoRequested = numDotsCyto(idx);
    if nucRequested > size(dotsNucleoid, 1)
        warning('Requested nucleoid dots (%d) exceed available nucleoid dots (%d). Using maximum available.', ...
            nucRequested, size(dotsNucleoid, 1));
        nucRequested = size(dotsNucleoid, 1);
    end
    if cytoRequested > size(dotsCytoplasm, 1)
        warning('Requested cytoplasmic dots (%d) exceed available cytoplasmic dots (%d). Using maximum available.', ...
            cytoRequested, size(dotsCytoplasm, 1));
        cytoRequested = size(dotsCytoplasm, 1);
    end

    % Randomly choose dots from nucleoid and cytoplasm.
    nucIndices = randperm(size(dotsNucleoid, 1), nucRequested);
    cytoIndices = randperm(size(dotsCytoplasm, 1), cytoRequested);
    coordForAccessibility{idx} = [dotsNucleoid(nucIndices, :); dotsCytoplasm(cytoIndices, :)];
end

% Generate the heat maps for each accessibility range setting.
cmpAccHeatMap = cell(length(accRange), 1);
for idx = 1:length(accRange)
    selectedDots = coordForAccessibility{idx};
    cmpAccHeatMap{idx} = Acc_htMap(selectedDots, scalingWidth, scalingLength, numPixY, aspectRatio, 'off');
end

end

%% ----- Supporting Functions -----

function [nucVolumeRatio, dotCoordinates] = Cal_nuc_vol(numDots, R, L, anchorModi)
% Cal_nuc_vol simulates dots within a cell and estimates the nucleoid volume ratio.
%
% Inputs:
%   numDots    - Total number of dots to simulate.
%   R          - Radius parameter for the cell.
%   L          - Length offset for the cell.
%   anchorModi - Matrix where each row defines an anchor point for the nucleoid boundary.
%
% Outputs:
%   nucVolumeRatio - The fraction of dots that fall inside the nucleoid.
%   dotCoordinates - A numDots-by-4 matrix of dots with columns [x, y, z, regionFlag].
%
% The regionFlag is 2 for nucleoid and 0 for cytoplasm.

% Define the cell boundary function.
% For a given y, cellBoundaryX returns the allowable x-value.
cellBoundaryX = @(y) ( (y >= -(R+L) & y < -L) .* sqrt(max(0, R^2 - (y+L).^2)) + ...
                       (y >= -L & y <= L)   .* R + ...
                       (y > L & y <= (L+R)) .* sqrt(max(0, R^2 - (y-L).^2)) );
                   
% Define the nucleoid (nucleus) boundary function using a piecewise linear function.
nucBoundaryX = @(y) compute_nucleus_boundary(anchorModi, y);

% Pre-allocate the dot coordinates matrix.
dotCoordinates = zeros(numDots, 4);
pointsGenerated = 0;
iteration = 0;
maxIterations = numDots * 100;  % Safety limit to prevent infinite loops

% Generate random dots until numDots points are accepted that fall inside the cell.
while pointsGenerated < numDots && iteration < maxIterations
    iteration = iteration + 1;
    % Random y between -(L+R) and L+R.
    yRand = -(L+R) + 2*(L+R)*rand(1);
    % Random x and z between -R and R.
    xRand = -R + 2*R*rand(1);
    zRand = -R + 2*R*rand(1);
    
    % Check if point is inside the cell
    if InRegion(xRand, yRand, zRand, cellBoundaryX)
        pointsGenerated = pointsGenerated + 1;
        dotCoordinates(pointsGenerated, 1:3) = [xRand, yRand, zRand];
    end
end

if pointsGenerated < numDots
    error('Failed to generate the required number of dots within the cell region. Check cell parameters.');
end

% Mark dots that fall within the nucleoid.
for i = 1:numDots
    if InRegion(dotCoordinates(i,1), dotCoordinates(i,2), dotCoordinates(i,3), nucBoundaryX)
        dotCoordinates(i,4) = 2;  % nucleoid flag
    else
        dotCoordinates(i,4) = 0;  % cytoplasmic flag
    end
end

% Compute the nucleoid volume ratio as the fraction of dots inside the nucleoid.
nucVolumeRatio = sum(dotCoordinates(:,4) ~= 0) / numDots;
end


function eq = compute_nucleus_boundary(anchor, y)
% compute_nucleus_boundary computes the piecewise linear boundary value 
% for the nucleoid at a given y.
%
% Inputs:
%   anchor - A matrix of anchor points. Each row defines an anchor point.
%
%   y      - The y-value(s) at which to compute the boundary value.
%
% Output:
%   eq     - The computed boundary x-value(s) corresponding to the input y 
%            according to the piecewise linear interpolation defined by 
%            the anchor points.

% Preallocate a matrix 'aa' to store the slope and intercept for each segment.
% 'aa(i,1)' will store the slope between anchor(i) and anchor(i+1)
% 'aa(i,2)' will store the intercept for that line segment.
aa = zeros(height(anchor)-1, 2);

% Loop over each consecutive pair of anchor points to compute slopes and intercepts.
for i = 1:height(anchor)-1
    dot1 = anchor(i, :);       
    dot2 = anchor(i+1, :);    
    % Compute the slope for the line segment connecting dot1 and dot2.
    aa(i, 1) = ((dot2(2) - dot1(2)) / (dot2(1) - dot1(1)));
    % Compute the intercept for the line segment using the determinant method.
    aa(i, 2) = (dot2(1) * dot1(2) - dot1(1) * dot2(2)) / (dot2(1) - dot1(1));
end
clear i dot1 dot2  % Clear loop variables to avoid namespace conflicts.

% Initialize the output variable.
eq = 0;    
% Compute the piecewise function value by summing over contributions from each segment.
for i = 1:height(anchor)-1
    dot1 = anchor(i, :);       
    dot2 = anchor(i+1, :);     
    eq = eq + ((y >= dot1(1)) & (y < dot2(1))) .* (aa(i, 1) .* y + aa(i, 2));
end
end

function isInside = InRegion(x, y, z, boundaryFunc)
% InRegion determines whether the point (x, y, z) lies within a region defined
% by a boundary function.
%
% Inputs:
%   x, y, z     - Scalars specifying a point’s coordinates.
%   boundaryFunc- Function handle that returns the boundary x-value when given y.
%
% Outputs:
%   isInside    - Logical true if the point is inside the region, false otherwise.
%
% The point is inside if the absolute value of its transformed x-coordinate is less
% than the absolute value given by the boundary function evaluated at y.

[transX, transY, ~] = CoordTrans(x, y, z);
isInside = abs(transX) < abs(boundaryFunc(transY));
end

function [transX, transY, angleXY] = CoordTrans(x, y, z)
% CoordTrans performs a coordinate transformation on a 3D point.
%
% Inputs:
%   x, y, z - The original coordinates.
%
% Outputs:
%   transX - The distance from central axis.
%   transY - The unchanged y-coordinate.
%   angleXY- The angle in the x-y plane computed with atan2.
%
transY = y;
radius = sqrt(x^2 + z^2);
if abs(x) < eps
    if z == 0
        transX = 0;
    else
        transX = radius * sign(z);
    end
else
    transX = radius * sign(x);
end
% angleXY = atan2(z, x);
angleXY = atan(z/x); % angle plane to x-y palne in radius; + --> clockwise
end

function heatMap = Acc_htMap(coordinates, scalingWidth, scalingLength, numPixY, aspectRatio, figOption)
% Acc_htMap generates a 2D heat map from dot coordinates.
%
% Inputs:
%   coordinates  - Matrix of dot coordinates (at least two columns for x and y).
%   scalingWidth - Scaling factor for the x coordinates.
%   scalingLength- Scaling factor for the y coordinates.
%   numPixY      - Number of pixels in the y-direction for histogram binning.
%   aspectRatio  - Aspect ratio for the heat map.
%   figOption    - 'on' to display the heat map plot; 'off' otherwise.
%
% Outputs:
%   heatMap      - The computed heat map (2D matrix).
%
% The function projects coordinates onto the x-y plane, scales them, symmetrizes the data,
% bins the data into a 2D histogram, and optionally plots the result.

% Project dot coordinates onto x-y plane with appropriate scaling.
projXY = coordinates(:, 1:2);
projXY(:,1) = projXY(:,1) / scalingWidth;
projXY(:,2) = projXY(:,2) / scalingLength;
projXY = projXY / 2;  % Additional normalization

% Create symmetric copies so that the heat map is symmetric about both axes.
mirrorX = projXY;
mirrorX(:,1) = -projXY(:,1);
symCoords = [projXY; mirrorX; -mirrorX; -projXY];

% Define the bin edges for histogram2.
xEdges = -0.5 : 1/(numPixY * aspectRatio) : 0.5;
yEdges = -0.5 : 1/numPixY : 0.5;

% Optionally plot the heat map.
if ~strcmpi(figOption, 'off')
    figure;
    histogram2(symCoords(:,2), symCoords(:,1), 'XBinEdges', xEdges, 'YBinEdges', yEdges, ...
        'Normalization','probability', 'DisplayStyle','tile', 'FaceColor','flat', 'ShowEmptyBins','on');
    colormap(figOption);
    colorbar;
    xlim([-0.5 0.5]);
    ylim([-0.5 0.5]);
    xticks([-0.5 -0.3 -0.1 0.1 0.3 0.5]);
    yticks([-0.5 -0.3 -0.1 0.1 0.3 0.5]);
    xticklabels({'0','0.2','0.4','0.6','0.8','1'});
    yticklabels({'0','0.2','0.4','0.6','0.8','1'});
    pbaspect([aspectRatio 1 1]);
    title('Localization Heat Map from Simulation');
end

% Compute (but do not plot) the heat map histogram.
heatMapCounts = histcounts2(symCoords(:,2), symCoords(:,1), 'XBinEdges', xEdges, 'YBinEdges', yEdges, ...
    'Normalization','probability');
heatMap = heatMapCounts.';
end