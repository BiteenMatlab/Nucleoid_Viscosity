<<<<<<< HEAD
function [pop_peri_ratio,pop_core_ratio] = Pop_Peri_Core_Ratio(num_dot,acc,R,L,anchor_spMap,...
    scalingWidth,scalingLength,num_pix_y,aspt_ratio,Peri)

    [~,dot_coord]=Cal_nuc_vol(num_dot,R,L,anchor_spMap);
    projXY = dot_coord(:,1:2); %clo1--> short; col2--> long
    
    projXY(:,1) = projXY(:,1) / scalingWidth;
    projXY(:,2) = projXY(:,2) / scalingLength;
    projXY = projXY / 2;  
    
    projXY = projXY + 0.5;
    projXY(:,1) = projXY(:,1)*num_pix_y;
    projXY(:,2) = projXY(:,2)*num_pix_y*aspt_ratio;
    
    peri1 = Peri{1};
    peri2 = Peri{2};
    
    peridot = inpolygon(projXY(:,1),projXY(:,2),peri1(:,2),peri1(:,1));
    coredot = inpolygon(projXY(:,1),projXY(:,2),peri2(:,2),peri2(:,1));
    
    buff_core = dot_coord(coredot==1,:);
    a = find(peridot==1);
    b = find(coredot==0);
    c = intersect(a,b);
    buff_peri = dot_coord(c,:);
    
    pop_peri_ratio = height(buff_peri(buff_peri(:,4)==2,:))/height(buff_peri)*acc;
    pop_core_ratio = height(buff_core(buff_core(:,4)==2,:))/height(buff_core)*acc;
end

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
=======
function [pop_peri_ratio,pop_core_ratio] = Pop_Peri_Core_Ratio(num_dot,acc,R,L,anchor_spMap,...
    scalingWidth,scalingLength,num_pix_y,aspt_ratio,Peri)

    [~,dot_coord]=Cal_nuc_vol(num_dot,R,L,anchor_spMap);
    projXY = dot_coord(:,1:2); %clo1--> short; col2--> long
    
    projXY(:,1) = projXY(:,1) / scalingWidth;
    projXY(:,2) = projXY(:,2) / scalingLength;
    projXY = projXY / 2;  
    
    projXY = projXY + 0.5;
    projXY(:,1) = projXY(:,1)*num_pix_y;
    projXY(:,2) = projXY(:,2)*num_pix_y*aspt_ratio;
    
    peri1 = Peri{1};
    peri2 = Peri{2};
    
    peridot = inpolygon(projXY(:,1),projXY(:,2),peri1(:,2),peri1(:,1));
    coredot = inpolygon(projXY(:,1),projXY(:,2),peri2(:,2),peri2(:,1));
    
    buff_core = dot_coord(coredot==1,:);
    a = find(peridot==1);
    b = find(coredot==0);
    c = intersect(a,b);
    buff_peri = dot_coord(c,:);
    
    pop_peri_ratio = height(buff_peri(buff_peri(:,4)==2,:))/height(buff_peri)*acc;
    pop_core_ratio = height(buff_core(buff_core(:,4)==2,:))/height(buff_core)*acc;
end

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
>>>>>>> 5e4cf7b (Initial commit)
end