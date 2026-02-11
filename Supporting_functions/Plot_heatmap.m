function [raw_heatmap, htMap] = Plot_heatmap(Store_track, Store_mask, Index, gap_y, aspt_ratio, fig, symmetrize, clr_map)
% Plot_heatmap Generates a 2D heat map (normalized localization 2D histogram) from
% cell track and mask data.
%
% [raw_heatmap, htMap] = Plot_heatmap(Store_track, Store_mask, Index, gap_y, aspt_ratio, fig, symmetrize, clr_map)
%
% The function first normalizes the localization coordinates of each cell using
% the corresponding cell shape (obtained from regionprops and feretProperties).
% It then optionally forces symmetry (by reflecting coordinates) and generates a
% 2D histogram (heat map). If requested, a figure is displayed.
%
% INPUTS:
% Store_track - Cell array. Each element contains a numeric matrix with track
% data for a ROI. It is assumed that column 5 holds cell IDs, and
% columns 2:3 contain localization coordinates.
% Store_mask - Cell array. Each element is a numeric matrix of the corresponding
% cell masks for a ROI.
% Index - Cell array. Each element is a vector of cell IDs (indices)
% corresponding to cells that pass the morphology filter.
% gap_y - (Optional) Scalar. Y‐axis bin width for the histogram.
% Default = 1/20.
% aspt_ratio - (Optional) Scalar. The aspect ratio used for normalization
% (x and y scaling). Default = 2.
% fig - (Optional) String, 'on' to display the heat map figure, or
% 'off' to skip plotting. Default = 'off'.
% symmetrize - (Optional) String, 'on' to force symmetry in localizations, or
% 'off' to leave as is. Default = 'on'.
% clr_map - (Optional) String specifying the colormap (e.g., 'parula').
% Default = 'parula'.
%
% OUTPUTS:
% raw_heatmap - A combined matrix (Nx2) of normalized localization coordinates
% from all ROIs.
% htMap - The 2D histogram (heat map matrix) produced by histcounts2.
%
% This function processes each ROI in parallel (using parfor) and, for each cell:
% 1. Extracts the single‐cell mask from the full-phase mask.
% 2. Extracts the corresponding track localizations.
% 3. Computes cell properties (e.g., via feretProperties) to obtain the cell
% center, dimensions, and orientation.
% 4. Converts the raw localization coordinates into normalized coordinates
% relative to cell dimensions.
% 5. Optionally removes outlier localizations.
% Finally, it combines these normalized coordinates, and, depending on the symmetrize
% flag, reflects the data to enforce symmetry before plotting the heat map.
%
% Author: Xiaofeng Dai
% Date: 05/28/2024

%% 1. Set Default Parameters and Validate Inputs
narginchk(3,8);
if nargin < 4, gap_y = 1/20; end
if nargin < 5, aspt_ratio = 2; end
if nargin < 6, fig = 'off'; end
if nargin < 7, symmetrize = 'on'; end
if nargin < 8, clr_map = 'parula'; end

if ~iscell(Store_track) || ~iscell(Store_mask) || ~iscell(Index)
    error('Store_track, Store_mask, and Index must be cell arrays.');
end

%% 2. Process Each ROI in Parallel
numFiles = numel(Store_track);
Result_all = cell(numFiles,1);

% Turn off warnings on parallel workers.
parfevalOnAll(gcp(), @warning, 0, 'off', 'all');

parfor i = 1:numFiles
    % Get current ROI data.
    tracks = Store_track{i};
    phaseMask = Store_mask{i};
    cellIDs = Index{i};  % Vector of cell IDs for this ROI.
    numCells = numel(cellIDs);
    
    % Preallocate cell arrays for each cell's mask and track localizations.
    cellMaskStore = cell(numCells,1);
    cellTrackStore = cell(numCells,1);
    
    if ~isempty(tracks)
        % Extract single‐cell masks.
        for j = 1:numCells
            currentID = cellIDs(j);
            singleCellMask = phaseMask;
            singleCellMask(singleCellMask ~= currentID) = 0;
            cellMaskStore{j} = singleCellMask;
        end
        % Extract corresponding localization coordinates.
        for j = 1:numCells
            currentID = cellIDs(j);
            locIndices = (tracks(:,5) == currentID);
            cellTrackStore{j} = tracks(locIndices, 2:3);
        end
    end
    
    % Process each cell within the ROI.
    resultSingleCell = cell(numCells,1);
    try
        for j = 1:numCells
            if ~isempty(cellTrackStore{j})
                % Get region properties of the single-cell mask.
                imgProps = regionprops('table', cellMaskStore{j}, 'PixelList');
                feretResults = feretProperties(imgProps);
                
                % Remove entries with empty PixelList.
                pixelLists = feretResults.PixelList;
                validIdx = false(numel(pixelLists),1);
                for k = 1:numel(pixelLists)
                    if ~isempty(pixelLists{k})
                        validIdx(k) = true;
                    end
                end
                feretResults = feretResults(validIdx,:);
                
                % Compute cell center using Minimum Area Bounding Box.
                boundingBox = cell2mat(feretResults.MinAreaBoundingBox);
                cellCenter = (boundingBox(1,:) + boundingBox(3,:)) / 2;
                
                % Extract cell dimension and orientation.
                xLength = feretResults.MaxFeretDiameter;
                if feretResults.MaxFeretDiameterOrientation > 0
                    poleAngleX = (180 - feretResults.MaxFeretDiameterOrientation) * pi/180;
                else
                    poleAngleX = - (pi * feretResults.MaxFeretDiameterOrientation) / 180;
                end
                yLength = feretResults.MinFeretDiameter;
                if feretResults.MinFeretDiameterOrientation > 0
                    poleAngleY = (180 - feretResults.MinFeretDiameterOrientation) * pi/180;
                else
                    poleAngleY = - (pi * feretResults.MinFeretDiameterOrientation) / 180;
                end
                
                % Process track localizations.
                trackData = cellTrackStore{j};
                numLocs = size(trackData,1);
                vectorData = zeros(numLocs,8);
                
                % Compute vector from cell center to each localization.
                vectorData(:,1) = trackData(:,2) - cellCenter(1);
                vectorData(:,2) = -trackData(:,1) + cellCenter(2);
                
                % Use a standard x-axis vector for angle calculation.
                baseVector = [1, 0];
                for k = 1:numLocs
                    vectorData(k,4) = anglevec(baseVector, [vectorData(k,1), vectorData(k,2)]);
                end
                
                % Compute Euclidean distance from the center.
                vectorData(:,3) = sqrt(vectorData(:,1).^2 + vectorData(:,2).^2);
                
                % Compute angular differences relative to pole angles.
                for k = 1:numLocs
                    thetaX = poleAngleX - vectorData(k,4);
                    if (thetaX <= -pi) && (thetaX > -2*pi)
                        vectorData(k,5) = mod(thetaX, pi);
                    else
                        vectorData(k,5) = thetaX;
                    end
                end
                for k = 1:numLocs
                    thetaY = poleAngleY - vectorData(k,4);
                    if (thetaY <= -pi) && (thetaY > -2*pi)
                        vectorData(k,6) = mod(thetaY, pi);
                    else
                        vectorData(k,6) = thetaY;
                    end
                end
                
                % Calculate normalized coordinates using the reference angle.
                angleRef = abs(poleAngleX - poleAngleY);
                vectorData(:,7) = (vectorData(:,3)/sin(angleRef)) .* sin(vectorData(:,6));
                vectorData(:,8) = (vectorData(:,3)/sin(angleRef)) .* sin(vectorData(:,5));
                vectorData(:,7) = vectorData(:,7) / xLength;
                vectorData(:,8) = vectorData(:,8) / yLength;
                
                % Store normalized coordinates.
                storeData = zeros(numLocs, 3);
                storeData(:,1:2) = vectorData(:,7:8);
                resultSingleCell{j} = storeData;
            end
        end
    catch ME
        warning('Error processing ROI %d: %s', i, ME.message);
    end
    
    % Remove outlier localizations for each cell.
    for j = 1:numCells
        if ~isempty(resultSingleCell{j})
            cellData = resultSingleCell{j};
            numPoints = size(cellData,1);
            outlierIdx1 = find(abs(cellData(:,1)) > 0.5);
            outlierIdx2 = find(abs(cellData(:,2)) > 0.5);
            totalOutliers = numel(outlierIdx1) + numel(outlierIdx2);
            if totalOutliers < 0.025 * numPoints
                % Remove outliers if they comprise less than 2.5% of the data.
                allOutliers = unique([outlierIdx1; outlierIdx2]);
                cellData(allOutliers,:) = [];
                resultSingleCell{j} = cellData;
            else
                % If too many points are outliers, discard this cell's data.
                resultSingleCell{j} = [];
            end
        end
    end
    
    % Concatenate results for the current ROI.
    resultFile = cat(1, resultSingleCell{:});
    Result_all{i} = resultFile;
end 

%% 3. Combine Data Across ROIs
raw_heatmap = cat(1, Result_all{:});

%% 4. Generate and (Optionally) Display Heat Map
gap_x = gap_y / aspt_ratio;
if strcmpi(fig, 'on')
    switch lower(symmetrize)
        case 'off'
            figure;
            histogram2(raw_heatmap(:,1), raw_heatmap(:,2), ...
                'XBinEdges', -0.5:gap_x:0.5, 'YBinEdges', -0.5:gap_y:0.5, ...
                'Normalization', 'probability', 'DisplayStyle', 'tile', ...
                'FaceColor', 'flat', 'ShowEmptyBins', 'on');
            colormap(gca, clr_map);
            colorbar;
            xlim([-0.5, 0.5]);
            xticks([-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]);
            ylim([-0.5, 0.5]);
            yticks([-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]);
            xticklabels({'0','0.2','0.4','0.6','0.8','1'});
            yticklabels({'0','0.2','0.4','0.6','0.8','1'});
            pbaspect([aspt_ratio, 1, 1]);
            
            testmap = histcounts2(raw_heatmap(:,1), raw_heatmap(:,2), ...
                'XBinEdges', -0.5:gap_x:0.5, 'YBinEdges', -0.5:gap_y:0.5, ...
                'Normalization', 'probability');
            htMap = testmap.';
            
        case 'on'
            % Force symmetry by reflecting the localizations.
            resultCat = raw_heatmap;
            for r = 1:size(resultCat,1)
                if resultCat(r,1) > 0
                    resultCat(r,1) = -resultCat(r,1);
                end
                if resultCat(r,2) > 0
                    resultCat(r,2) = -resultCat(r,2);
                end
            end
            % Create mirrored copies.
            resultCat3 = [resultCat(:,1), -resultCat(:,2)];
            resultCat4 = [-resultCat(:,1), -resultCat(:,2)];
            resultCat5 = [-resultCat(:,1), resultCat(:,2)];
            heatmap_sym = [resultCat(:,1:2); resultCat3; resultCat4; resultCat5];
            
            figure;
            histogram2(heatmap_sym(:,1), heatmap_sym(:,2), ...
                'XBinEdges', -0.5:gap_x:0.5, 'YBinEdges', -0.5:gap_y:0.5, ...
                'Normalization', 'probability', 'DisplayStyle', 'tile', ...
                'FaceColor', 'flat', 'ShowEmptyBins', 'on');
            colormap(gca, clr_map);
            colorbar;
            xlim([-0.5, 0.5]);
            xticks([-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]);
            ylim([-0.5, 0.5]);
            yticks([-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]);
            xticklabels({'0','0.2','0.4','0.6','0.8','1'});
            yticklabels({'0','0.2','0.4','0.6','0.8','1'});
            pbaspect([aspt_ratio, 1, 1]);
            
            testmap = histcounts2(heatmap_sym(:,1), heatmap_sym(:,2), ...
                'XBinEdges', -0.5:gap_x:0.5, 'YBinEdges', -0.5:gap_y:0.5, ...
                'Normalization', 'probability');
            htMap = testmap.';
            
        otherwise
            error('Invalid option for symmetrize. Use ''on'' or ''off''.');
    end
else
    % If no figure should be displayed, compute htMap without plotting.
    switch lower(symmetrize)
        case 'off'
            testmap = histcounts2(raw_heatmap(:,1), raw_heatmap(:,2), ...
                'XBinEdges', -0.5:gap_x:0.5, 'YBinEdges', -0.5:gap_y:0.5, ...
                'Normalization', 'probability');
            htMap = testmap.';
        case 'on'
            resultCat = raw_heatmap;
            for r = 1:size(resultCat,1)
                if resultCat(r,1) > 0
                    resultCat(r,1) = -resultCat(r,1);
                end
                if resultCat(r,2) > 0
                    resultCat(r,2) = -resultCat(r,2);
                end
            end
            resultCat3 = [resultCat(:,1), -resultCat(:,2)];
            resultCat4 = [-resultCat(:,1), -resultCat(:,2)];
            resultCat5 = [-resultCat(:,1), resultCat(:,2)];
            heatmap_sym = [resultCat(:,1:2); resultCat3; resultCat4; resultCat5];
            
            testmap = histcounts2(heatmap_sym(:,1), heatmap_sym(:,2), ...
                'XBinEdges', -0.5:gap_x:0.5, 'YBinEdges', -0.5:gap_y:0.5, ...
                'Normalization', 'probability');
            htMap = testmap.';
        otherwise
            error('Invalid option for symmetrize. Use ''on'' or ''off''.');
    end
end
end