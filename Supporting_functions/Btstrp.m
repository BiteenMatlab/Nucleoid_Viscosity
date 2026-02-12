<<<<<<< HEAD
function bootstrp_list = Btstrp(cellIndices, numSample, numIter, randomSeed)
% Btstrp performs bootstrap sampling on cell indices.
%
% bootstrp_list = Btstrp(cellIndices, numSample, numIter) returns a cell
% array with bootstrap samples. The input cellIndices is a cell array,
% where each cell contains a numeric vector of cell IDs (e.g., those passing
% a morphology filter) from a given ROI. This function first builds a combined
% list with two columns:
% Column 1: ROI index
% Column 2: Cell ID (from cellIndices)
%
% Then, for each of the numIter iterations, it randomly samples (with replacement)
% numSample rows from the combined list.
%
% bootstrp_list = Btstrp(cellIndices, numSample, numIter, randomSeed) allows you to
% specify a random seed (randomSeed) for reproducibility. If randomSeed is provided,
% the function sets the RNG to that seed for the duration of bootstrap sampling and
% then restores the previous RNG state. If not provided or empty, the function uses the
% current global RNG state.
%
% INPUTS:
% cellIndices - Cell array where each element is a numeric vector of cell IDs
% for each ROI.
% numSample - Positive scalar specifying how many cells to sample in each bootstrap iteration.
% numIter - Positive scalar specifying the number of bootstrap iterations.
% randomSeed - (Optional) Scalar specifying the seed to use for RNG.
%
% OUTPUT:
% bootstrp_list - A cell array (numIter x 1), where each element is a numeric
% matrix of size [numSample x 2] containing [ROI_index, cell_ID].
%
% Example:
% % Using a fixed seed for reproducibility:
% bootList = Btstrp(cellIndices, 100, 500, 1);
%
% Author: Xiaofeng Dai
% Date: 04/11/2025

%% 1. Input Validation
if ~iscell(cellIndices)
    error('cellIndices must be a cell array.');
end
if nargin < 3
    error('Three inputs (cellIndices, numSample, numIter) are required.');
end
if ~isscalar(numSample) || numSample <= 0
    error('numSample must be a positive scalar.');
end
if ~isscalar(numIter) || numIter <= 0
    error('numIter must be a positive scalar.');
end

%% 2. Optionally Set the Random Seed
if nargin >= 4 && ~isempty(randomSeed)
    oldState = rng;  % Save the current RNG state
    rng(randomSeed); % Set the RNG to randomSeed
end

%% 3. Build Combined List of ROI and Cell IDs
% Each row: [ROI_index, cell_ID]
numROIs = numel(cellIndices);
listID = cell(numROIs, 1);
for roiIdx = 1:numROIs
    currentCellIDs = cellIndices{roiIdx};
    if ~isempty(currentCellIDs)
        currentCellIDs = currentCellIDs(:);  % Ensure a column vector
        listID{roiIdx} = [repmat(roiIdx, numel(currentCellIDs), 1), currentCellIDs];
    end
end
listID = vertcat(listID{:});
numCells = size(listID, 1);
% check if numSample is greater than total number of cells
if numSample > numCells
    error('number of samples is larger than the total number of cells %d',numel(vertcat(cellIndices{:})))
end
%% 4. Bootstrap Sampling
bootstrp_list = cell(numIter, 1);
for iter = 1:numIter
    sample_indices = randsample(numCells, numSample, true);
    bootstrp_list{iter} = listID(sample_indices, :);
end

%% 5. Restore the Original RNG State if a Seed Was Used
if exist('oldState', 'var')
    rng(oldState);
end
=======
function bootstrp_list = Btstrp(cellIndices, numSample, numIter, randomSeed)
% Btstrp performs bootstrap sampling on cell indices.
%
% bootstrp_list = Btstrp(cellIndices, numSample, numIter) returns a cell
% array with bootstrap samples. The input cellIndices is a cell array,
% where each cell contains a numeric vector of cell IDs (e.g., those passing
% a morphology filter) from a given ROI. This function first builds a combined
% list with two columns:
% Column 1: ROI index
% Column 2: Cell ID (from cellIndices)
%
% Then, for each of the numIter iterations, it randomly samples (with replacement)
% numSample rows from the combined list.
%
% bootstrp_list = Btstrp(cellIndices, numSample, numIter, randomSeed) allows you to
% specify a random seed (randomSeed) for reproducibility. If randomSeed is provided,
% the function sets the RNG to that seed for the duration of bootstrap sampling and
% then restores the previous RNG state. If not provided or empty, the function uses the
% current global RNG state.
%
% INPUTS:
% cellIndices - Cell array where each element is a numeric vector of cell IDs
% for each ROI.
% numSample - Positive scalar specifying how many cells to sample in each bootstrap iteration.
% numIter - Positive scalar specifying the number of bootstrap iterations.
% randomSeed - (Optional) Scalar specifying the seed to use for RNG.
%
% OUTPUT:
% bootstrp_list - A cell array (numIter x 1), where each element is a numeric
% matrix of size [numSample x 2] containing [ROI_index, cell_ID].
%
% Example:
% % Using a fixed seed for reproducibility:
% bootList = Btstrp(cellIndices, 100, 500, 1);
%
% Author: Xiaofeng Dai
% Date: 04/11/2025

%% 1. Input Validation
if ~iscell(cellIndices)
    error('cellIndices must be a cell array.');
end
if nargin < 3
    error('Three inputs (cellIndices, numSample, numIter) are required.');
end
if ~isscalar(numSample) || numSample <= 0
    error('numSample must be a positive scalar.');
end
if ~isscalar(numIter) || numIter <= 0
    error('numIter must be a positive scalar.');
end

%% 2. Optionally Set the Random Seed
if nargin >= 4 && ~isempty(randomSeed)
    oldState = rng;  % Save the current RNG state
    rng(randomSeed); % Set the RNG to randomSeed
end

%% 3. Build Combined List of ROI and Cell IDs
% Each row: [ROI_index, cell_ID]
numROIs = numel(cellIndices);
listID = cell(numROIs, 1);
for roiIdx = 1:numROIs
    currentCellIDs = cellIndices{roiIdx};
    if ~isempty(currentCellIDs)
        currentCellIDs = currentCellIDs(:);  % Ensure a column vector
        listID{roiIdx} = [repmat(roiIdx, numel(currentCellIDs), 1), currentCellIDs];
    end
end
listID = vertcat(listID{:});
numCells = size(listID, 1);
% check if numSample is greater than total number of cells
if numSample > numCells
    error('number of samples is larger than the total number of cells %d',numel(vertcat(cellIndices{:})))
end
%% 4. Bootstrap Sampling
bootstrp_list = cell(numIter, 1);
for iter = 1:numIter
    sample_indices = randsample(numCells, numSample, true);
    bootstrp_list{iter} = listID(sample_indices, :);
end

%% 5. Restore the Original RNG State if a Seed Was Used
if exist('oldState', 'var')
    rng(oldState);
end
>>>>>>> 5e4cf7b (Initial commit)
end