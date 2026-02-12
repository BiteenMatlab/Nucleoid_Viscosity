function [coordInDomain, scalingCoord, hFig] = FindDomainCoord( ...
                spMap, macroOutline, nucOutline, ...
                numPixY, asptRatio, figPlot)
%--------------------------------------------------------------------------
% FindDomainCoord
% Identifies pixel coordinates that lie inside both
%   (a) the macrodomain region, and
%   (b) the nucleoid region.
%
% INPUTS
%   spMap        : sptail displacement map
%   macroOutline : P×2 array – [x y] outline points of macrodomain region
%   nucOutline   : Q×2 array – [x y] outline points of nucleoid
%   numPixY      : scalar    – #pixels along the y-axis
%   asptRatio    : scalar    – x:y pixel size ratio of the map
%   figPlot      : 'on'/'off' (optional, default 'off')
%
% OUTPUTS
%   coordInDomain : N×2 array – all valid (x,y) pixel indices
%   domainPoly    : cell      – coordinates from each sub region after scaling
%   hFig          : figure handle (empty if figPlot = 'off')
%--------------------------------------------------------------------------

%% --------------------------- 1. Defaults --------------------------------
if nargin < 6 || isempty(figPlot)
    figPlot = 'off';
end

%% --------------------------- 2. Validation ------------------------------
assert(ismatrix(spMap)           , '`spMap` must be a 2-D numeric matrix.');
assert(size(macroOutline,2) == 2 , '`macroOutline` must be P×2.');
assert(size(nucOutline  ,2) == 2 , '`nucOutline`   must be Q×2.');
assert(numPixY   > 0 && asptRatio > 0, 'numPixY and asptRatio must be > 0.');

%% --------------------------- 3. Cluster macro outline -------------------
% DBSCAN groups contiguous macro points into sub-regions (noise is idx = -1)
idx = dbscan(macroOutline, 0.1, 5);
clusterIDs = unique(idx(idx > 0));        % discard the unassigned

if isempty(clusterIDs)
    warning('No macrodomain sub-regions identified by DBSCAN.');
    coordInDomain = zeros(0,2);
    scalingCoord    = {};
    hFig          = [];
    return
end

nRegions   = numel(clusterIDs);
scalingCoord = cell(nRegions,1);            % keep polygon vertices per region
coordCell  = cell(nRegions,1);            % temp storage for coordinates

%% --------------------------- 4. Generate full-grid coords ---------------
[gridX, gridY] = meshgrid(1:width(spMap), 1:height(spMap));
mapCoord      = [gridX(:), gridY(:)];     % (#pixels) × 2

%% --------------------------- 5. Pre-compute nucleoid mask ---------------
inNuc = inpolygon(mapCoord(:,1), mapCoord(:,2), ...
                  nucOutline(:,1), nucOutline(:,2));

%% --------------------------- 6. Loop over sub-regions -------------------
for k = 1:nRegions
    %--- 6.1  Current sub-region vertices ---------------------------------
    subReg  = macroOutline(idx == clusterIDs(k), :);

    % boundary() returns vertices of a tight boundary
    bndIdx  = boundary(subReg(:,1), subReg(:,2), 0.6);
    poly    = subReg(bndIdx, :);          % 2-column array

    %--- 6.2  Convert (normalised) outline → pixel coordinates ------------
    poly(:,1) = (poly(:,1) - 0.5) * numPixY           + numPixY/2 + 0.5;          % y-axis
    poly(:,2) = (poly(:,2) - 0.5) * numPixY*asptRatio + numPixY*asptRatio/2 + 0.5; % x-axis
    scalingCoord{k} = poly;                 % store for output / plotting

    %--- 6.3  Points inside current domain *and* nucleoid -----------------
    inDom = inpolygon(mapCoord(:,1), mapCoord(:,2), ...
                      poly(:,2), poly(:,1));           % note swapped columns
    coordCell{k} = mapCoord(inDom & inNuc, :);
end

%% --------------------------- 7. Concatenate results ---------------------
coordInDomain = vertcat(coordCell{:});

%% --------------------------- 8. Plot if requested -----------------------
if strcmpi(figPlot,'on')
    oldState = warning('off','MATLAB:polyshape:repairedBySimplify');
    hFig = figure('Name','Domain & nucleoid pixels');
    % Nucleoid outline
    plot(polyshape(nucOutline(1:end-1,1), nucOutline(1:end-1,2)), ...
         'FaceColor', '#D95319', 'EdgeAlpha', 0.05, 'EdgeColor', '#D95319');
    hold on
    % All selected pixels
    scatter(coordInDomain(:,1), coordInDomain(:,2), 3, 'b', 'filled');

    % Each domain polygon (scaled)
    for k = 1:nRegions
        poly = scalingCoord{k};
        plot(polyshape(poly(1:end-1,2), poly(1:end-1,1)), ...
             'FaceColor', '#30E8EB', 'EdgeAlpha', 0.05, 'EdgeColor', '#30E8EB');
    end
    axis equal
    % xlabel('x (pixel index)'); ylabel('y (pixel index)');
    xlim([0 numPixY*asptRatio])
    ylim([0 numPixY])
    title('Pixels inside macrodomain ∩ nucleoid');
    warning(oldState);
else
    hFig = [];
end
end

% function pts = stripDuplicates(pts)
% % Remove consecutive duplicate vertices and any closing vertex that equals
% % the first.  Keeps the original order (important for polygons).
%     if isequal(pts(1,:), pts(end,:))
%         pts(end,:) = [];
%     end
%     keep = [true; any(diff(pts),2)];   % diff == 0 → duplicate row
%     pts  = pts(keep,:);
% end