function [cellMorph, morphByMovie] = CalMorph(pixelSize, maskStore, figOption, numBins)
% CalMorph calculates cell morphological characteristics from segmented cell masks.
%
% This function processes a cell mask of each region of interest (ROI) to compute
% the length and width of individual cells using Feret diameter analysis (function 
% feretProperties in Folder 'Feret_analysis' by Steve Eddins @ Mathwork). It returns 
% both a concatenated matrix of morphology metrics and a cell array with per-ROI results.
%
% INPUTS:
% pixelSize - (Scalar, numeric) Pixel size in nanometers (nm/pixel). Must be > 0.
% maskStore - (Cell array) Each element is a 2D numeric array (ROI) where
% nonzero integer labels denote individual cells.
% figOption - (Optional, string) 'on' to display histogram plots; 'off' (default) otherwise.
% numBins - (Optional, integer) Number of bins for histograms (default 50 when figOption='on').
%
% OUTPUTS:
% cellMorph - (Matrix) Combined cell morphology measures for all ROIs. Each row corresponds to one cell,
% and the columns represent:
% 1) Cell ID (in its ROI)
% 2) Long Axis length (µm)
% 3) Short Axis length (µm)
% 4) Aspect Ratio (Long/Short)
% 5) Cell Area (µm^2)
% morphByMovie - (Cell array) Each cell element is the morphology matrix (as above) for one ROI.
%
% NOTE: The feretProperties function (which must be on your MATLAB path) is assumed to take as input
% the output of regionprops with 'PixelList' and return a structure with
% fields 'MaxFeretDiameter' and 'MinFeretDiameter'.
%
% EXAMPLE:
% [cellMorph, morphByMovie] = CalMorph(49, Store_mask, 'on', 50);
%
% Author: Xiaofeng Dai
% Date: 05/24/2024

%% 1. Input Validation and Defaults
narginchk(2, 4);

% Validate pixelSize
if ~isnumeric(pixelSize) || ~isscalar(pixelSize) || pixelSize <= 0
    error('pixelSize must be a positive numeric scalar.');
end

% Validate maskStore as cell array
if ~iscell(maskStore)
    error('maskStore must be a cell array of segmented mask matrices.');
end

% Default for figOption if not provided
if nargin < 3 || isempty(figOption)
    figOption = 'off';
end

% Default for numBins if not provided (only needed if figure plotting is on)
if nargin < 4 || isempty(numBins)
    if strcmpi(figOption, 'on')
        numBins = 50;
    end
end

%% 2. Preallocate Result Container for ROIs
numROIs = numel(maskStore);
morphResultsCell = cell(numROIs, 1);

%% 3. Process Each ROI (Parallelized)
% Each ROI is processed to extract individual cell morphology.
parfor roiIdx = 1:numROIs           
    currentMask = maskStore{roiIdx};
    % Verify that the current ROI mask is numeric and nonempty.
    if isempty(currentMask) || ~isnumeric(currentMask)
        morphResultsCell{roiIdx} = [];
        continue;
    end
    
    % Determine the number of cells based on the maximum label in the mask.
    numCells = max(currentMask, [], 'all');
    % Preallocate matrix with 5 columns:
    % Column1: Cell ID, Column2: Long Axis (µm), Column3: Short Axis (µm),
    % Column4: Aspect Ratio, Column5: Area (µm^2)
    cellMorphology = zeros(numCells, 5);
    cellMorphology(:, 1) = (1:numCells)';     % Assign cell IDs
    
    roiSize = size(currentMask);
    
    % Process each cell in the ROI
    for cellIdx = 1:numCells
        % Create a logical (binary) mask for the current cell.
        binaryCellMask = (currentMask == cellIdx);
        if ~any(binaryCellMask(:))
            % Skip if there are no pixels for this label (should rarely happen)
            continue;
        end
        
        % Get pixel coordinates for the cell.
        [rowIdx, colIdx] = find(binaryCellMask);
        
        % Define a bounding box around the cell with a 10-pixel margin.
        rowMin = max(min(rowIdx) - 10, 1);
        rowMax = min(max(rowIdx) + 10, roiSize(1));
        colMin = max(min(colIdx) - 10, 1);
        colMax = min(max(colIdx) + 10, roiSize(2));
        
        % Crop the binary mask to the bounding box.
        croppedCell = binaryCellMask(rowMin:rowMax, colMin:colMax);
        
        % Calculate region properties (only PixelList is needed for feret analysis).
        cellPropsTable = regionprops('table', croppedCell, 'PixelList');

        % if cell mask is consisted of discontiguous regions
        if height(cellPropsTable.PixelList)>1
            % % focus on the largest contiguous region
            % heights = cellfun(@(x) height(x), cellPropsTable.PixelList);
            % [~,longest_id] = max(heights);
            % id_range = 1:1:height(cellPropsTable.PixelList);
            % id_range(longest_id) = [];
            % cellPropsTable(id_range,:) = [];

            % ignore this cell
            warning('discontiguous regions in cell %d in ROI %d.', cellIdx, roiIdx);
            continue;
        end

        % Compute Feret diameters via the feretProperties function.
        feretResults = feretProperties(cellPropsTable);
        
        % Extract Feret diameters.
        longAxis = feretResults.MaxFeretDiameter;
        shortAxis = feretResults.MinFeretDiameter;
        aspectRatio = longAxis / shortAxis;

        % Convert from pixels to micrometers.
        % (Assuming pixelSize is in nm so that: pixel in µm = pixelSize/1000)
        cellMorphology(cellIdx, 2) = longAxis * pixelSize / 1000;   % Long axis in µm
        cellMorphology(cellIdx, 3) = shortAxis * pixelSize / 1000;  % Short axis in µm
        cellMorphology(cellIdx, 4) = aspectRatio;                   % Aspect ratio (dimensionless)
        % Compute area. (Number of pixels * (pixelSize^2 in nm^2)) converted to µm^2.
        cellMorphology(cellIdx, 5) = nnz(binaryCellMask) * (pixelSize^2) / 1e6;
    end
    morphResultsCell{roiIdx} = cellMorphology;
end

%% 4. Consolidate and Filter Results
% Store per-ROI results.
morphByMovie = morphResultsCell;
% Concatenate all ROI results vertically.
allMorphResults = vertcat(morphResultsCell{:});
% Keep only cells with a nonzero long axis measurement.
cellMorph = allMorphResults(allMorphResults(:,2) ~= 0, :);

%% 5. Plot Histogram Figures (if enabled)
if strcmpi(figOption, 'on')
    figure;
    figure('Units','pixels','Position',[100 100 1200 800]);
    % Histogram of Cell Long Axis
    subplot(2,2,1);
    meanLong = mean(cellMorph(:,2));
    stdLong  = std(cellMorph(:,2));
    numLong = numel(cellMorph(:,2));
    histTitle = sprintf('%.2f ± %.2f μm (N = %d)', meanLong, stdLong, numLong);
    histogram(cellMorph(:,2), numBins, 'Normalization', 'probability');
    title({'Length of Long Axis', histTitle}, 'FontSize', 15);
    xlabel('Cell Length (μm)');
    ylabel('Probability');
    
    % Histogram of Cell Short Axis
    subplot(2,2,2);
    meanShort = mean(cellMorph(:,3));
    stdShort  = std(cellMorph(:,3));
    numShort = numel(cellMorph(:,3));
    histTitle = sprintf('%.2f ± %.2f μm (N = %d)', meanShort, stdShort, numShort);
    histogram(cellMorph(:,3), numBins, 'Normalization', 'probability');
    title({'Length of Short Axis', histTitle}, 'FontSize', 15);
    xlabel('Cell Width (μm)');
    ylabel('Probability');
    
    % Histogram of Aspect Ratio
    subplot(2,2,3);
    meanAspect = mean(cellMorph(:,4));
    stdAspect  = std(cellMorph(:,4));
    numAspect  = numel(cellMorph(:,4));
    histTitle = sprintf('%.2f ± %.2f (N = %d)', meanAspect, stdAspect, numAspect);
    histogram(cellMorph(:,4), numBins, 'Normalization', 'probability');
    title({'Aspect Ratio', histTitle}, 'FontSize', 15);
    xlabel('Aspect Ratio');
    ylabel('Probability');
    
    % Histogram of Cell Area
    subplot(2,2,4);
    meanArea = mean(cellMorph(:,5));
    stdArea  = std(cellMorph(:,5));
    numArea  = numel(cellMorph(:,5));
    histTitle = sprintf('%.2f ± %.2f μm² (N = %d)', meanArea, stdArea, numArea);
    histogram(cellMorph(:,5), numBins, 'Normalization', 'probability');
    title({'Cell Area', histTitle}, 'FontSize', 15);
    xlabel('Cell Area (μm²)');
    ylabel('Probability');
end
end