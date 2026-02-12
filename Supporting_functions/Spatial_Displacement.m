function [raw_spatialDisp, Map, Count_Map] = Spatial_Displacement(Store_track, Store_mask, Index, pixel_size, min_steps, gap_y, aspt_ratio, thres_ratio, fig, symmetrize, clr_map)
% Spatial_Displacement Creates a spatial displacement map from localization data.
%
% [raw_spatialDisp, Map] = Spatial_Displacement(Store_track, Store_mask, Index, pixel_size, ...
% min_steps, gap_y, aspt_ratio, thres_ratio, fig, symmetrize, clr_map)
%
% This function first normalizes the coordinates of each localization based on the
% shape of the corresponding cells. Then, for each cell, it calculates the local step
% sizes between consecutive localizations (if a minimum number of steps is reached for a trajectory).
% Finally, the function bins these normalized coordinates into a grid and computes
% the average step size (in nm) at each grid location.
%
% Optional inputs (if not provided, defaults are used):
% pixel_size - Pixel size (default 49, unit nm).
% min_steps - Minimum number of steps for analysis per trajectory (default 3).
% gap_y - Bin width on y-axis for the grid (default 1/20).
% aspt_ratio - Aspect ratio used for spatial map (default 2).
% thres_ratio - Threshold ratio for filtering low count bins (default 0.1).
% fig - 'on' to display figures, 'off' otherwise (default 'off').
% symmetrize - 'on' forces symmetry of the data, 'off' keeps original orientation (default 'on').
% clr_map - Colormap for plotting (default 'jet').
%
% OUTPUTS:
% raw_spatialDisp - A matrix [N x 3] where columns 1 and 2 are the normalized
% x and y coordinates and column 3 is the step size (in physical units).
% Map - The final grid (heat map) of average displacement, where bins
% with low counts have been set to NaN.
% Count_Map - Grid map with number of steps from each specific bin.
%
% Author: Xiaofeng Dai
% Date: 06/01/2024

%% 1. Set Default Parameters for Optional Inputs
if nargin < 4, pixel_size  = 49;      end
if nargin < 5, min_steps   = 3;        end
if nargin < 6, gap_y       = 1/20;     end
if nargin < 7, aspt_ratio  = 2;        end
if nargin < 8, thres_ratio = 0.1;      end
if nargin < 9, fig         = 'off';    end
if nargin < 10, symmetrize = 'on';     end
if nargin < 11, clr_map    = 'jet';    end

% Validate that the three main inputs are cell arrays.
if ~iscell(Store_track) || ~iscell(Store_mask) || ~iscell(Index)
    error('Store_track, Store_mask, and Index must be cell arrays.');
end

%% 2. Process Each ROI in Parallel
numFiles = numel(Store_track);
SSS_result = cell(numFiles, 1);

% Turn off warnings on parallel workers
parfevalOnAll(gcp(), @warning, 0, 'off', 'all');

parfor i = 1:numFiles        
    % Retrieve data for ROI i.
    tracks   = Store_track{i};
    phaseMsk = Store_mask{i};
    cellIDs  = Index{i};
    numCells = numel(cellIDs);
    
    % Preallocate cell arrays for single-cell masks and track data.
    cellMaskStore = cell(numCells, 1);
    cellTrackStore = cell(numCells, 1);
    
    if ~isempty(tracks)
        % Extract the mask for each cell.
        for j = 1:numCells
            currID = cellIDs(j);
            singleCell = (phaseMsk(:) == currID);
            cellMaskStore{j} = reshape(singleCell, size(phaseMsk));
        end
        
        % Extract corresponding localizations (assuming column 5 holds cell ID)
        for j = 1:numCells
            currID = cellIDs(j);
            % Use logical indexing instead of looping to extract rows.
            cellTrackStore{j} = tracks(tracks(:,5)==currID, :);
        end
    end
    
    % For each cell, calculate the step sizes in its track, if available.
    StepSize_Store = cell(numCells, 1);
    for j = 1:numCells
        if ~isempty(cellTrackStore{j})
            % Extract unique track numbers (assumed stored in column 4)
            trackNums = unique(cellTrackStore{j}(:,4));
            stepsPerTrack = cell(numel(trackNums), 1);
            for k = 1:numel(trackNums)
                % Use direct indexing for extraction.
                idx = cellTrackStore{j}(:,4) == trackNums(k);
                % Use helper function ArrayExtract to get the rows.
                trackSegment = ArrayExtract(cellTrackStore{j}, find(idx));
                if size(trackSegment, 1) > min_steps
                    % Calculate step sizes for this track.
                    sss_tmp = SpatialStepSize(trackSegment);
                    stepsPerTrack{k} = sss_tmp;
                end
            end
            % Concatenate step size data for the cell.
            cellSteps = cat(1, stepsPerTrack{:});
            % Remove rows that are all zeros.
            cellSteps(all(cellSteps==0,2), :) = [];
            StepSize_Store{j} = cellSteps;
        end
    end
    
    % Compute normalized displacement for each cell.
    Result_Singlecell = cell(numCells, 1);
    for j = 1:numCells
        if ~isempty(StepSize_Store{j})
            % Compute cell shape properties from the single-cell mask.
            imgProps    = regionprops('table', cellMaskStore{j}, 'PixelList');
            feretResult = feretProperties(imgProps);
            
            % For discontinuous segmented results, select the region with max amount of pixel as cell region.
            if numel(feretResult.PixelList) > 1
                counts = zeros(numel(feretResult.PixelList), 1);
                for k = 1:numel(feretResult.PixelList)
                    if ~isempty(feretResult.PixelList{k})
                        counts(k) = size(feretResult.PixelList{k}, 1);
                    end
                end
                % Select rows with the maximum pixel count.
                feretResult = feretResult(counts == max(counts), :);
            end
            
            % Compute cell center from the minimum-area bounding box.
            bbox   = cell2mat(feretResult.MinAreaBoundingBox);
            center = (bbox(1,:) + bbox(3,:)) * 0.5;
            Xlength = feretResult.MaxFeretDiameter;
            % Compute PoleAngleX (convert degrees to radians)
            if feretResult.MaxFeretDiameterOrientation > 0
                PoleAngleX = (180 - feretResult.MaxFeretDiameterOrientation)*pi/180;
            else
                PoleAngleX = -(pi * feretResult.MaxFeretDiameterOrientation)/180;
            end
            Ylength = feretResult.MinFeretDiameter;
            if feretResult.MinFeretDiameterOrientation > 0
                PoleAngleY = (180 - feretResult.MinFeretDiameterOrientation)*pi/180;
            else
                PoleAngleY = -(pi * feretResult.MinFeretDiameterOrientation)/180;
            end
            
            tmpSteps = StepSize_Store{j};
            numLocs  = size(tmpSteps, 1);
            vectData = zeros(numLocs, 8);
            % Column 1-2: Vector from cell center to localization (normalized later)
            vectData(:,1) = tmpSteps(:,2) - center(1);
            vectData(:,2) = -tmpSteps(:,1) + center(2);
            baseVec = [1, 0];
            % Column 4: Angle (in radians) of the vector relative to x-axis.
            for k = 1:numLocs
                vectData(k,4) = anglevec(baseVec, [vectData(k,1), vectData(k,2)]);
            end
            % Column 3: Euclidean distance.
            vectData(:,3) = sqrt(vectData(:,1).^2 + vectData(:,2).^2);
            % Column 5: Difference from PoleAngleX.
            for k = 1:numLocs
                thetaX = PoleAngleX - vectData(k,4);
                if (thetaX <= -pi) && (thetaX > -2*pi)
                    vectData(k,5) = mod(thetaX, pi);
                else
                    vectData(k,5) = thetaX;
                end
            end
            % Column 6: Difference from PoleAngleY.
            for k = 1:numLocs
                thetaY = PoleAngleY - vectData(k,4);
                if (thetaY <= -pi) && (thetaY > -2*pi)
                    vectData(k,6) = mod(thetaY, pi);
                else
                    vectData(k,6) = thetaY;
                end
            end
            % Calculate normalized coordinates.
            angleRef = abs(PoleAngleX - PoleAngleY);
            vectData(:,7) = (vectData(:,3)/sin(angleRef)) .* sin(vectData(:,6));
            vectData(:,8) = (vectData(:,3)/sin(angleRef)) .* sin(vectData(:,5));
            % Normalize by cell dimensions.
            vectData(:,7) = vectData(:,7) / Xlength;
            vectData(:,8) = vectData(:,8) / Ylength;
            
            % Store normalized coordinates and step size.
            storeData = zeros(numLocs, 3);
            storeData(:,1:2) = vectData(:,7:8);
            storeData(:,3) = tmpSteps(:,3);
            Result_Singlecell{j} = storeData;
        end
    end
    
    % Outlier removal for each cell:
    for j = 1:numCells
        if ~isempty(Result_Singlecell{j})
            cellData = Result_Singlecell{j};
            numPoints = size(cellData,1);
            outlierIdx1 = find(abs(cellData(:,1)) > 0.5);
            outlierIdx2 = find(abs(cellData(:,2)) > 0.5);
            totalOut = numel(outlierIdx1) + numel(outlierIdx2);
            if isempty(outlierIdx1) && isempty(outlierIdx2)
                % Nothing to remove.
            elseif totalOut < 0.025 * numPoints
                combinedOut = unique([outlierIdx1; outlierIdx2]);
                cellData(combinedOut,:) = [];
                Result_Singlecell{j} = cellData;
            else
                % If too many outliers, discard this cell's data.
                Result_Singlecell{j} = [];
            end
        end
    end
    % Concatenate results for this ROI.
    ROI_result = cat(1, Result_Singlecell{:});
    SSS_result{i} = ROI_result;
end 

%% 3. Combine Data and Scale Step Size
Result_cat = cat(1, SSS_result{:});
% Convert step size from pixel units to physical dimensions.
Result_cat(:,3) = Result_cat(:,3) * pixel_size;
raw_spatialDisp = Result_cat;

%% 4. Generate Spatial Displacement Map
% Define grid dimensions.
num_pix_y = 1 / gap_y;
num_pix_x = aspt_ratio * num_pix_y;

if strcmpi(fig, 'on')
    switch lower(symmetrize)
        case 'off'
            % Compute grid using helper function PlotSSS.
            [DisZ, CntZ] = PlotSSS(Result_cat, num_pix_y, aspt_ratio);
            % Determine threshold from top 10% of count bins.
            countsVec = reshape(CntZ, [], 1);
            topCounts = maxk(countsVec, round(0.1 * length(countsVec)));
            threshold = thres_ratio * mean(topCounts);
            % Set bins below threshold to NaN.
            DisZ(CntZ < threshold) = NaN;
            Count_Map = CntZ;
            Map = DisZ;
            figure;
            imagesc(DisZ);
            axis equal;
            colormap(gca, clr_map);
            colorbar;
            xlim([0.5, num_pix_x+0.5]);
            xticks(linspace(0.5, num_pix_x+0.5, 6));
            ylim([0.5, num_pix_y+0.5]);
            yticks(linspace(0.5, num_pix_y+0.5, 6));
            xticklabels({'0','0.2','0.4','0.6','0.8','1'});
            yticklabels({'1','0.8','0.6','0.4','0.2','0'});
            pbaspect([aspt_ratio, 1, 1]);
            
        case 'on'
            [DisZ, CntZ] = PlotSSS(Result_cat, num_pix_y, aspt_ratio * 2);
            SumDisZ = DisZ .* CntZ;
            Map_Shape = size(DisZ);
            % Reflect quadrants to force symmetry.
            Region1 = SumDisZ(1:Map_Shape(1)/2, 1:Map_Shape(2)/2);  % upper-left
            Count1  = CntZ(1:Map_Shape(1)/2, 1:Map_Shape(2)/2);
            Region2 = SumDisZ((Map_Shape(1)/2+1):Map_Shape(1), 1:Map_Shape(2)/2); % bottom-left
            Count2  = CntZ((Map_Shape(1)/2+1):Map_Shape(1), 1:Map_Shape(2)/2);
            Region2 = flipud(Region2); Count2 = flipud(Count2);
            Region3 = SumDisZ(1:Map_Shape(1)/2, (Map_Shape(2)/2+1):Map_Shape(2)); % upper-right
            Count3  = CntZ(1:Map_Shape(1)/2, (Map_Shape(2)/2+1):Map_Shape(2));
            Region3 = fliplr(Region3); Count3 = fliplr(Count3);
            Region4 = SumDisZ((Map_Shape(1)/2+1):Map_Shape(1), (Map_Shape(2)/2+1):Map_Shape(2)); % bottom-right
            Count4  = CntZ((Map_Shape(1)/2+1):Map_Shape(1), (Map_Shape(2)/2+1):Map_Shape(2));
            Region4 = rot90(Region4,2); Count4 = rot90(Count4,2);
            Regions = Region1 + Region2 + Region3 + Region4;
            Counts  = Count1 + Count2 + Count3 + Count4;
            Map1 = Regions ./ Counts;
            Map2 = flipud(Map1);
            Map3 = fliplr(Map1);
            Map4 = rot90(Map1,2);
            Map_combined = [Map1, Map3; Map2, Map4];
            Count_Map = [Counts, fliplr(Counts); flipud(Counts), rot90(Counts,2)];
            
            % Average pair-wise columns.
            buff_map = zeros(size(Map_combined,1), size(Map_combined,2)/2);
            buff_count = zeros(size(Map_combined,1), size(Map_combined,2)/2);
            for cc = 1:(size(Map_combined,2)/2)
                col1 = Map_combined(:, 2*cc - 1);
                col2 = Map_combined(:, 2*cc);
                cnt1 = Count_Map(:, 2*cc - 1);
                cnt2 = Count_Map(:, 2*cc);
                buff_map(:, cc) = (col1 .* cnt1 + col2 .* cnt2) ./ (cnt1 + cnt2);
                buff_count(:, cc) = (cnt1 + cnt2);
            end
            Count_Map = buff_count;
            Map_combined = buff_map;

            % Set bins with low counts to NaN.
            countsVec = reshape(Count_Map, [], 1);
            topCounts = maxk(countsVec, round(0.1 * length(countsVec)));
            threshold = thres_ratio * mean(topCounts);
            Map_combined(Count_Map < threshold) = NaN;
            Map = Map_combined;
            figure;
            imagesc(Map);
            axis equal;
            colormap(gca, clr_map);
            colorbar;
            xlim([0.5, num_pix_x+0.5]);
            xticks(linspace(0.5, num_pix_x+0.5, 6));
            ylim([0.5, num_pix_y+0.5]);
            yticks(linspace(0.5, num_pix_y+0.5, 6));
            xticklabels({'0','0.2','0.4','0.6','0.8','1'});
            yticklabels({'1','0.8','0.6','0.4','0.2','0'});
            pbaspect([aspt_ratio, 1, 1]);
    end
else
    % If no figure is requested, compute Map without plotting.
    switch lower(symmetrize)
        case 'off'
            [DisZ, CntZ] = PlotSSS(Result_cat, num_pix_y, aspt_ratio);
            DisZ(CntZ < (thres_ratio * mean(maxk(reshape(CntZ,[],1), round(0.1*numel(CntZ))))) ) = NaN;
            Map = DisZ;
            Count_Map = CntZ;
        case 'on'
            [DisZ, CntZ] = PlotSSS(Result_cat, num_pix_y, aspt_ratio*2);
            SumDisZ = DisZ .* CntZ;
            Map_Shape = size(DisZ);
            Region1 = SumDisZ(1:(Map_Shape(1)/2), 1:(Map_Shape(2)/2));
            Count1  = CntZ(1:(Map_Shape(1)/2), 1:(Map_Shape(2)/2));
            Region2 = SumDisZ((Map_Shape(1)/2+1):Map_Shape(1), 1:(Map_Shape(2)/2));
            Count2  = CntZ((Map_Shape(1)/2+1):Map_Shape(1), 1:(Map_Shape(2)/2));
            Region2 = flipud(Region2); Count2 = flipud(Count2);
            Region3 = SumDisZ(1:(Map_Shape(1)/2), (Map_Shape(2)/2+1):Map_Shape(2));
            Count3  = CntZ(1:(Map_Shape(1)/2), (Map_Shape(2)/2+1):Map_Shape(2));
            Region3 = fliplr(Region3); Count3 = fliplr(Count3);
            Region4 = SumDisZ((Map_Shape(1)/2+1):Map_Shape(1), (Map_Shape(2)/2+1):Map_Shape(2));
            Count4  = CntZ((Map_Shape(1)/2+1):Map_Shape(1), (Map_Shape(2)/2+1):Map_Shape(2));
            Region4 = rot90(Region4,2); Count4 = rot90(Count4,2);
            Regions = Region1 + Region2 + Region3 + Region4;
            Counts  = Count1 + Count2 + Count3 + Count4;
            Map1 = Regions ./ Counts;
            Map2 = flipud(Map1);
            Map3 = fliplr(Map1);
            Map4 = rot90(Map1,2);
            Map_combined = [Map1, Map3; Map2, Map4];
            Count_Map = [Counts, fliplr(Counts); flipud(Counts), rot90(Counts,2)];
            
            buff_map = zeros(size(Map_combined,1), size(Map_combined,2)/2);
            buff_count = zeros(size(Map_combined,1), size(Map_combined,2)/2);
            for cc = 1:(size(Map_combined,2)/2)
                col1 = Map_combined(:, 2*cc-1);
                col2 = Map_combined(:, 2*cc);
                cnt1 = Count_Map(:, 2*cc-1);
                cnt2 = Count_Map(:, 2*cc);
                buff_map(:, cc) = (col1.*cnt1 + col2.*cnt2) ./ (cnt1+cnt2);
                buff_count(:, cc) = (cnt1+cnt2);
            end
            Count_Map = buff_count;
            Map_combined = buff_map;
            countsVec = reshape(Count_Map, [], 1);
            topCounts = maxk(countsVec, round(0.1*numel(Count_Map)));
            threshold = thres_ratio * mean(topCounts);
            Map_combined(Count_Map < threshold) = NaN;
            Map = Map_combined;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Function: SpatialStepSize
% Calculates the step size between consecutive localizations. Only considers
% steps when the first column difference equals -1 (e.g., consecutive frames).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SSS = SpatialStepSize(trackinfo)
numRows = size(trackinfo, 1);
SSS = zeros(numRows - 1, 3);
for i = 1:(numRows - 1)
% Check if the frame difference is exactly -1.
if (trackinfo(i,1) - trackinfo(i+1,1)) == -1
% Compute the midpoint coordinates and the Euclidean distance.
midRow = 0.5 * (trackinfo(i,2) + trackinfo(i+1,2));
midCol = 0.5 * (trackinfo(i,3) + trackinfo(i+1,3));
step = sqrt((trackinfo(i,2)-trackinfo(i+1,2))^2 + (trackinfo(i,3)-trackinfo(i+1,3))^2);
SSS(i, :) = [midRow, midCol, step];
end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Function: ArrayExtract
% Extracts rows from a matrix based on the provided indices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function buffer = ArrayExtract(inputMatrix, indices)
buffer = inputMatrix(indices, :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Function: PlotSSS
% Bins the normalized step size data into a grid and computes the average step size.
% NormSSS is an [N x 3] matrix: columns 1 and 2 are normalized x and y coordinates,
% column 3 is the step size.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DisZ, CntZ] = PlotSSS(NormSSS, n, aspect_ratio)
xscale = 1 / (aspect_ratio * n);
yscale = 1 / n;
DisZ = zeros(n, aspect_ratio * n);
CntZ = zeros(n, aspect_ratio * n);
numPoints = size(NormSSS,1);
for i = 1:numPoints
% Determine grid indices for each localization.
xIndex = floor((NormSSS(i,1) + 0.5) / xscale) + 1;
yIndex = floor((NormSSS(i,2) + 0.5) / yscale) + 1;
% Accumulate step size and count.
DisZ(yIndex, xIndex) = DisZ(yIndex, xIndex) + NormSSS(i,3);
CntZ(yIndex, xIndex) = CntZ(yIndex, xIndex) + 1;
end
% Compute average displacement per bin; avoid division by zero.
DisZ(CntZ > 0) = DisZ(CntZ > 0) ./ CntZ(CntZ > 0);
DisZ(CntZ == 0) = NaN;
end