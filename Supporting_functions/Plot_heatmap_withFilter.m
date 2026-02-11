function [raw_heatmap, htMap, mean_loc_prec] = Plot_heatmap_withFilter(Store_track, Store_mask, Store_fits, Index, gap_y, aspt_ratio, uncertainty, fig, symmetrize, clr_map)
% Plot_heatmap_withFilter Computes a localization heat map with an additional
% localization precision filter.
%
% [raw_heatmap, htMap, mean_loc_prec] = Plot_heatmap_withFilter(Store_track, Store_mask, Store_fits, Index, gap_y, aspt_ratio, uncertainty, fig, symmetrize, clr_map)
%
% This function processes tracking data over multiple ROIs. In each ROI, the localization
% confidence intervals (widthcCI) from the Store_fits structure are used to filter out poor
% localizations based on the specified uncertainty threshold. The remaining localizations are
% then used to compute a normalized localization heat map (raw_heatmap) and a binned heat map (htMap).
% Finally, the average localization precision (in pixels) is computed from the fitted data.
%
% Inputs:
% Store_track - Cell array containing tracking data matrices for each ROI.
% Store_mask - Cell array of corresponding cell mask matrices.
% Store_fits - Cell array of fit structures for each ROI. Each structure must contain
% The 3 cells above are results from SMALLLABS tracking analysis.
% fields 'goodfit' (logical vector) and 'widthcCI' (numeric vector of CIs).
% Index - Cell array; each element is a vector of cell IDs for that ROI.
% gap_y - (Optional) Scalar bin width for the y-axis (default = 1/20).
% aspt_ratio - (Optional) Scalar indicating cell aspect ratio (default = 2).
% uncertainty - (Optional) Numeric threshold for widthcCI filtering (default = 2.5).
% fig - (Optional) 'on' to display figures; 'off' otherwise (default if not provided see below).
% symmetrize - (Optional) 'on' to force symmetry in the heat map; 'off' otherwise.
% clr_map - (Optional) Colormap for the plot (default = 'parula').
%
% Outputs:
% raw_heatmap - An [N x 3] matrix of normalized localization information.
% htMap - 2D histogram (heat map) matrix computed from raw_heatmap.
% mean_loc_prec - Mean localization precision (in pixels) computed from the fit data.
%
% Example:
% [rawHM, htMap, mlp] = Plot_heatmap_withFilter(Store_track, Store_mask, Store_fits, Index, 1/20, 2, 100, 'on', 'on', 'parula');
%
% Author: Xiaofeng Dai
% Date: 05/06/2025

%% 1. Set Default Parameters
if nargin < 10 || isempty(clr_map)
    clr_map = 'parula';
end
if nargin < 9 || isempty(symmetrize)
    symmetrize = 'on';
end
if nargin < 8 || isempty(fig)
    fig = 'on';
end
if nargin < 7 || isempty(uncertainty)
    uncertainty = 2.5;
end
if nargin < 6 || isempty(aspt_ratio)
    aspt_ratio = 2;
end
if nargin < 5 || isempty(gap_y)
    gap_y = 1/20;
end

%% 2. Process Each ROI in Parallel to Compute Normalized Localizations
numFiles = numel(Store_track);
Result_all = cell(numFiles,1);

% Disable warnings on all parallel workers.
parfevalOnAll(gcp(), @warning, 0, 'off', 'all');

parfor i = 1:numFiles
    % Get current ROI tracking data.
    tracks = Store_track{i};

    % ** Apply Localization Precision Filter **
    % Retrieve the "goodfit" flag and "widthcCI" from the fit structure.
    goodfit = Store_fits{i}.goodfit;
    widthcCI = Store_fits{i}.widthcCI;
    % col 6 in tracks is molecule indices.
    numTracks = size(tracks, 1);
    validTrackFlag = zeros(numTracks, 1);
    for iter = 1:numTracks
        molID = tracks(iter,6);
        % Accept the localization if goodfit is true and the CI is below uncertainty.
        if goodfit(molID) == 1 && widthcCI(molID) < uncertainty
            validTrackFlag(iter) = 1;
        end
    end
    tracks = tracks(validTrackFlag == 1, :);

    % Retrieve the cell mask and cell IDs for the ROI.
    phaseMask = Store_mask{i};
    cellIDs = Index{i};
    numCells = numel(cellIDs);

    % Preallocate cell arrays to store each cell's mask and localization coordinates.
    cellMaskStore = cell(numCells, 1);
    cellTrackStore = cell(numCells, 1);

    if ~isempty(tracks)
        % Get the individual cell mask for each cell.
        for j = 1:numCells
            currentID = cellIDs(j);
            singleCellMask = phaseMask;
            singleCellMask(singleCellMask ~= currentID) = 0;
            cellMaskStore{j} = singleCellMask;
        end
        % Extract the corresponding localization coordinates.
        for j = 1:numCells
            currentID = cellIDs(j);
            locInd = (tracks(:,5) == currentID);
            cellTrackStore{j} = tracks(locInd, 2:3);
        end
    end

    % Process each cell's localizations.
    resultSingleCell = cell(numCells, 1);
    try
        for j = 1:numCells
            if ~isempty(cellTrackStore{j})
                % Compute region properties using the cell mask.
                imgProps = regionprops('table', cellMaskStore{j}, 'PixelList');
                feretResults = feretProperties(imgProps);
                
                % Remove entries with empty PixelList.
                pixelLists = feretResults.PixelList;
                validIdx = false(numel(pixelLists), 1);
                for k = 1:numel(pixelLists)
                    if ~isempty(pixelLists{k})
                        validIdx(k) = true;
                    end
                end
                feretResults = feretResults(validIdx, :);
                
                % Compute cell center using the Minimum Area Bounding Box method.
                boundingBox = cell2mat(feretResults.MinAreaBoundingBox);
                cellCenter = (boundingBox(1,:) + boundingBox(3,:)) / 2;
                
                % Extract cell dimensions (MaxFeretDiameter & MinFeretDiameter) and orientations.
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

                % Process the tracking localizations.
                trackData = cellTrackStore{j};
                numLocs = size(trackData, 1);
                vectorData = zeros(numLocs, 8);
                % Compute vectors from the cell center.
                vectorData(:,1) = trackData(:,2) - cellCenter(1);
                vectorData(:,2) = -trackData(:,1) + cellCenter(2);
                baseVector = [1, 0];
                for k = 1:numLocs
                    vectorData(k,4) = anglevec(baseVector, [vectorData(k,1), vectorData(k,2)]);
                end
                vectorData(:,3) = sqrt(vectorData(:,1).^2 + vectorData(:,2).^2);

                % Calculate angular differences relative to pole angles.
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

                % Compute normalized coordinates based on the reference angle.
                angleRef = abs(poleAngleX - poleAngleY);
                vectorData(:,7) = (vectorData(:,3) / sin(angleRef)) .* sin(vectorData(:,6));
                vectorData(:,8) = (vectorData(:,3) / sin(angleRef)) .* sin(vectorData(:,5));
                vectorData(:,7) = vectorData(:,7) / xLength;
                vectorData(:,8) = vectorData(:,8) / yLength;

                storeData = zeros(numLocs, 3);
                storeData(:,1:2) = vectorData(:,7:8);
                resultSingleCell{j} = storeData;
            end
        end
    catch ME
        warning('Error processing ROI %d: %s', i, ME.message);
    end

    % Outlier removal for each cell.
    for j = 1:numCells
        if ~isempty(resultSingleCell{j})
            cellData = resultSingleCell{j};
            numPoints = size(cellData, 1);
            outIdx1 = find(abs(cellData(:,1)) > 0.5);
            outIdx2 = find(abs(cellData(:,2)) > 0.5);
            totalOut = numel(outIdx1) + numel(outIdx2);
            if totalOut < 0.025 * numPoints
                cellData(unique([outIdx1; outIdx2]), :) = [];
                resultSingleCell{j} = cellData;
            else
                resultSingleCell{j} = [];
            end
        end
    end
    
    % Concatenate results for this ROI.
    roiResult = cat(1, resultSingleCell{:});
    Result_all{i} = roiResult;
end

%% 3. Combine Data Across ROIs
raw_heatmap = cat(1, Result_all{:});

%% 4. Generate and (Optionally) Display the Heat Map
gap_x = gap_y / aspt_ratio;
if strcmpi(fig, 'on')
    switch lower(symmetrize)
        case 'off'
            figure;
            histogram2(raw_heatmap(:,1), raw_heatmap(:,2), ...
                'XBinEdges', -0.5:gap_x:0.5, 'YBinEdges', -0.5:gap_y:0.5, ...
                'Normalization', 'probability', 'DisplayStyle', 'tile', ...
                'FaceColor', 'flat', 'ShowEmptyBins', 'on','EdgeColor','none');
            colormap(gca, clr_map);
            colorbar;
            xlim([-0.5, 0.5]);
            xticks([-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]);
            ylim([-0.5, 0.5]);
            yticks([-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]);
            xticklabels({'0','0.2','0.4','0.6','0.8','1'});
            yticklabels({'0','0.2','0.4','0.6','0.8','1'});
            pbaspect([aspt_ratio, 1, 1]);
            % title('Localization Heatmap')
            
            testmap = histcounts2(raw_heatmap(:,1), raw_heatmap(:,2), ...
                'XBinEdges', -0.5:gap_x:0.5, 'YBinEdges', -0.5:gap_y:0.5, ...
                'Normalization', 'probability');
            htMap = testmap.';
            
        case 'on'
            % Force symmetry by reflecting the localization data.
            resultCat = raw_heatmap;
            for r = 1:size(resultCat,1)
                if resultCat(r,1) > 0, resultCat(r,1) = -resultCat(r,1); end
                if resultCat(r,2) > 0, resultCat(r,2) = -resultCat(r,2); end
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
                'FaceColor', 'flat', 'ShowEmptyBins', 'on','EdgeColor','none');
            colormap(gca, clr_map);
            colorbar;
            xlim([-0.5, 0.5]);
            xticks([-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]);
            ylim([-0.5, 0.5]);
            yticks([-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]);
            xticklabels({'0','0.2','0.4','0.6','0.8','1'});
            yticklabels({'0','0.2','0.4','0.6','0.8','1'});
            pbaspect([aspt_ratio, 1, 1]);
            % title('Localization Heatmap')
            
            testmap = histcounts2(heatmap_sym(:,1), heatmap_sym(:,2), ...
                'XBinEdges', -0.5:gap_x:0.5, 'YBinEdges', -0.5:gap_y:0.5, ...
                'Normalization', 'probability');
            htMap = testmap.';
            
        otherwise
            error('Invalid option for symmetrize. Use ''on'' or ''off''.');
    end
else
    % Compute htMap without plotting.
    switch lower(symmetrize)
        case 'off'
            testmap = histcounts2(raw_heatmap(:,1), raw_heatmap(:,2), ...
                'XBinEdges', -0.5:gap_x:0.5, 'YBinEdges', -0.5:gap_y:0.5, ...
                'Normalization', 'probability');
            htMap = testmap.';
        case 'on'
            resultCat = raw_heatmap;
            for r = 1:size(resultCat,1)
                if resultCat(r,1) > 0, resultCat(r,1) = -resultCat(r,1); end
                if resultCat(r,2) > 0, resultCat(r,2) = -resultCat(r,2); end
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

%% 5. Compute Mean Localization Precision (Diagnostic)
% Extract the localization precision values from all fit structures.
allFits = cat(1, Store_fits{:});
fitFields = struct2cell(allFits);
locPrecField = fitFields(16, :);
locPrecAll = cat(1, locPrecField{:});
% Remove outliers that exceed the uncertainty threshold.
locPrecFiltered = locPrecAll(locPrecAll <= uncertainty);
% Remove NaNs.
locPrecFiltered(isnan(locPrecFiltered)) = [];
mean_loc_prec = mean(locPrecFiltered);
% fprintf('Mean localization precision (pixels): %.3f\n', mean_loc_prec);
end
