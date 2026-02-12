<<<<<<< HEAD
function Metric_optz = htMap_similarity(refHtMap, cmpAccHtMap, accRange, lb, rb, figPlot)
% htMap_similiraty computes similarity metrics between a reference heat
% map and a set of heat maps over a range of accessibility values.
%
%   Metric_optz = htMap_similiraty(refHtMap, cmpAccHtMap, accRange, lb, rb, figPlot)
%
% Inputs:
%   refHtMap     - Reference heat map matrix.
%   cmpAccHtMap  - Cell array of comparison heat map matrices.
%   accRange     - A vector of accessibility range factors (one per cell in cmpAccHtMap).
%   lb           - Left bound index (column index) for the region of interest.
%   rb           - Right bound index (column index) for the region of interest.
%   figPlot      - 'on' to display plots; 'off' to suppress plotting.
%
% Output:
%   Metric_optz - A matrix where each row corresponds to one value in accRange.
%                 The columns are:
%                      Column 1: L2 norm distance (lower is better)
%                      Column 2: Pearson correlation coefficient (higher is better)
%                      Column 3: SSIM (higher is better)
%                      
%
% Example:
%   metrics = htMap_similiraty(refHtMap, cmpAccHtMap, 0:0.1:1, 10, 30, 'on');

%% ----- Defensive Checks -----
% Check that the reference heat map is numeric.
if ~isnumeric(refHtMap)
    error('The reference heat map (refHtMap) must be numeric.');
end

% Check that cmpAccHtMap is a cell array.
if ~iscell(cmpAccHtMap)
    error('cmpAccHtMap must be a cell array of heat map matrices.');
end

% Check that the length of accRange matches the number of cells in cmpAccHtMap.
if numel(accRange) ~= numel(cmpAccHtMap)
    error('The length of accRange must equal the number of cell entries in cmpAccHtMap.');
end

% Ensure that lb and rb are valid indices for the heat map columns.
[numRows, numCols] = size(refHtMap);
if lb < 1 || rb > numCols || lb > rb
    error('Bounds lb and rb must be valid column indices within the reference heat map.');
end

% Validate the plotting flag input.
if ~ischar(figPlot) && ~isstring(figPlot)
    error('figPlot must be either ''on'' or ''off''.');
end

% Define a small constant to avoid log(0) issues.
epsVal = 1e-12;

%% ----- Extract and Normalize the Reference Heat Map Region  -----
% Extract the region defined by columns lb to rb.
refSegment = refHtMap(:, lb:rb);

% Check if the region sum is nonzero.
refSum = sum(refSegment, 'all');
if refSum == 0
    error('The sum of the selected region in the reference heat map is 0.');
end

% Normalize the extracted region to form a probability distribution.
refSegment = refSegment / refSum;

%% ----- Initialize Metric Arrays -----
numAcc = numel(accRange);
Metric_distance     = zeros(numAcc, 1);  % Squared L2 norm difference.
Metric_pcc          = zeros(numAcc, 1);  % Pearson correlation coefficient.
Metric_SSIM         = zeros(numAcc, 1);  % SSIM.
% Metric_JSD = zeros(numAcc, 1);  % Jensen-Shannon divergence.

%% ----- Loop through each accessibility value and compare heat maps -----
for idx = 1:numAcc
    % Get the current comparison heat map from the cell array.
    currentHtMap = cmpAccHtMap{idx};
    
    % Check that the current heat map is numeric and has enough columns.
    if ~isnumeric(currentHtMap)
        error('Each entry in cmpAccHtMap must be numeric.');
    end
    if ~isequal(size(currentHtMap) , size(refHtMap))
        error('A heat map in cmpAccHtMap does not match dimension of refHtMap.');
    end
    
    % Extract the region [lb:rb] from the comparison heat map.
    cmpSegment = currentHtMap(:, lb:rb);
    
    % Normalize the comparison segment.
    cmpSum = sum(cmpSegment, 'all');
    if cmpSum == 0
        error('The sum of the selected region in a comparison heat map is 0.');
    end
    cmpSegment = cmpSegment / cmpSum;
    
    % Compute L2 norm distance.
    Metric_distance(idx) = sum((refSegment - cmpSegment) .^ 2, 'all');
    
    % Compute Pearson Correlation Coefficient.
    % Reshape the matrices into vectors.
    refVec = refSegment(:);
    cmpVec = cmpSegment(:);
    pccMatrix = corrcoef(refSegment, cmpSegment);
    Metric_pcc(idx) = pccMatrix(1, 2);
    
    % Compute SSIM.
    Metric_SSIM(idx) = weightedSSIM(refSegment, cmpSegment);
    
    % % Compute log of JSD.
    % Metric_JSD(idx) = log(weightedJSD(refSegment, cmpSegment));
end

%% ----- Optional Plotting of the Metrics -----
if strcmpi(figPlot, 'on')
    % Determine tick spacing based on accessibility range span.
    if max(accRange) - min(accRange) > 0.5
        gap = 0.1;
    else
        gap = 0.02;
    end
    
    % Create a new figure and plot each metric.
    figure;
    figure('Units','pixels','Position',[100 100 1000 1500]);
    % Plot L2 norm distance.
    subplot(3, 1, 1)
    plot(accRange, Metric_distance, 'LineWidth', 2);
    title('L2-norm Distance (minimize)', 'FontSize', 18);
    xticks(min(accRange):gap:max(accRange));
    set(gca, 'LineWidth', 1.5);
    
    % Plot Pearson Correlation Coefficient.
    subplot(3, 1, 2)
    plot(accRange, Metric_pcc, 'LineWidth', 2);
    title('Pearson Correlation Coefficient (maximize)', 'FontSize', 18);
    xticks(min(accRange):gap:max(accRange));
    set(gca, 'LineWidth', 1.5);
    
    % Plot SSIM.
    subplot(3, 1, 3)
    plot(accRange, Metric_SSIM, 'LineWidth', 2);
    title('SSIM (maximize)', 'FontSize', 18);
    xticks(min(accRange):gap:max(accRange));
    set(gca, 'LineWidth', 1.5);
    
    % % Plot JS Divergence.
    % subplot(4, 1, 4)
    % plot(accRange, Metric_JSD, 'LineWidth', 2);
    % title('Log(JS Divergence) (minimize)', 'FontSize', 18);
    % xticks(min(accRange):gap:max(accRange));
    % set(gca, 'LineWidth', 1.5);
end

%% ----- Combine Metrics into Output Matrix -----
% Each row corresponds to an element in accRange with columns:
% [L2 distance, PCC, Cosine similarity]
Metric_optz = cat(2, Metric_distance, Metric_pcc, Metric_SSIM);

=======
function Metric_optz = htMap_similarity(refHtMap, cmpAccHtMap, accRange, lb, rb, figPlot)
% htMap_similiraty computes similarity metrics between a reference heat
% map and a set of heat maps over a range of accessibility values.
%
%   Metric_optz = htMap_similiraty(refHtMap, cmpAccHtMap, accRange, lb, rb, figPlot)
%
% Inputs:
%   refHtMap     - Reference heat map matrix.
%   cmpAccHtMap  - Cell array of comparison heat map matrices.
%   accRange     - A vector of accessibility range factors (one per cell in cmpAccHtMap).
%   lb           - Left bound index (column index) for the region of interest.
%   rb           - Right bound index (column index) for the region of interest.
%   figPlot      - 'on' to display plots; 'off' to suppress plotting.
%
% Output:
%   Metric_optz - A matrix where each row corresponds to one value in accRange.
%                 The columns are:
%                      Column 1: L2 norm distance (lower is better)
%                      Column 2: Pearson correlation coefficient (higher is better)
%                      Column 3: SSIM (higher is better)
%                      
%
% Example:
%   metrics = htMap_similiraty(refHtMap, cmpAccHtMap, 0:0.1:1, 10, 30, 'on');

%% ----- Defensive Checks -----
% Check that the reference heat map is numeric.
if ~isnumeric(refHtMap)
    error('The reference heat map (refHtMap) must be numeric.');
end

% Check that cmpAccHtMap is a cell array.
if ~iscell(cmpAccHtMap)
    error('cmpAccHtMap must be a cell array of heat map matrices.');
end

% Check that the length of accRange matches the number of cells in cmpAccHtMap.
if numel(accRange) ~= numel(cmpAccHtMap)
    error('The length of accRange must equal the number of cell entries in cmpAccHtMap.');
end

% Ensure that lb and rb are valid indices for the heat map columns.
[numRows, numCols] = size(refHtMap);
if lb < 1 || rb > numCols || lb > rb
    error('Bounds lb and rb must be valid column indices within the reference heat map.');
end

% Validate the plotting flag input.
if ~ischar(figPlot) && ~isstring(figPlot)
    error('figPlot must be either ''on'' or ''off''.');
end

% Define a small constant to avoid log(0) issues.
epsVal = 1e-12;

%% ----- Extract and Normalize the Reference Heat Map Region  -----
% Extract the region defined by columns lb to rb.
refSegment = refHtMap(:, lb:rb);

% Check if the region sum is nonzero.
refSum = sum(refSegment, 'all');
if refSum == 0
    error('The sum of the selected region in the reference heat map is 0.');
end

% Normalize the extracted region to form a probability distribution.
refSegment = refSegment / refSum;

%% ----- Initialize Metric Arrays -----
numAcc = numel(accRange);
Metric_distance     = zeros(numAcc, 1);  % Squared L2 norm difference.
Metric_pcc          = zeros(numAcc, 1);  % Pearson correlation coefficient.
Metric_SSIM         = zeros(numAcc, 1);  % SSIM.
% Metric_JSD = zeros(numAcc, 1);  % Jensen-Shannon divergence.

%% ----- Loop through each accessibility value and compare heat maps -----
for idx = 1:numAcc
    % Get the current comparison heat map from the cell array.
    currentHtMap = cmpAccHtMap{idx};
    
    % Check that the current heat map is numeric and has enough columns.
    if ~isnumeric(currentHtMap)
        error('Each entry in cmpAccHtMap must be numeric.');
    end
    if ~isequal(size(currentHtMap) , size(refHtMap))
        error('A heat map in cmpAccHtMap does not match dimension of refHtMap.');
    end
    
    % Extract the region [lb:rb] from the comparison heat map.
    cmpSegment = currentHtMap(:, lb:rb);
    
    % Normalize the comparison segment.
    cmpSum = sum(cmpSegment, 'all');
    if cmpSum == 0
        error('The sum of the selected region in a comparison heat map is 0.');
    end
    cmpSegment = cmpSegment / cmpSum;
    
    % Compute L2 norm distance.
    Metric_distance(idx) = sum((refSegment - cmpSegment) .^ 2, 'all');
    
    % Compute Pearson Correlation Coefficient.
    % Reshape the matrices into vectors.
    refVec = refSegment(:);
    cmpVec = cmpSegment(:);
    pccMatrix = corrcoef(refSegment, cmpSegment);
    Metric_pcc(idx) = pccMatrix(1, 2);
    
    % Compute SSIM.
    Metric_SSIM(idx) = weightedSSIM(refSegment, cmpSegment);
    
    % % Compute log of JSD.
    % Metric_JSD(idx) = log(weightedJSD(refSegment, cmpSegment));
end

%% ----- Optional Plotting of the Metrics -----
if strcmpi(figPlot, 'on')
    % Determine tick spacing based on accessibility range span.
    if max(accRange) - min(accRange) > 0.5
        gap = 0.1;
    else
        gap = 0.02;
    end
    
    % Create a new figure and plot each metric.
    figure;
    figure('Units','pixels','Position',[100 100 1000 1500]);
    % Plot L2 norm distance.
    subplot(3, 1, 1)
    plot(accRange, Metric_distance, 'LineWidth', 2);
    title('L2-norm Distance (minimize)', 'FontSize', 18);
    xticks(min(accRange):gap:max(accRange));
    set(gca, 'LineWidth', 1.5);
    
    % Plot Pearson Correlation Coefficient.
    subplot(3, 1, 2)
    plot(accRange, Metric_pcc, 'LineWidth', 2);
    title('Pearson Correlation Coefficient (maximize)', 'FontSize', 18);
    xticks(min(accRange):gap:max(accRange));
    set(gca, 'LineWidth', 1.5);
    
    % Plot SSIM.
    subplot(3, 1, 3)
    plot(accRange, Metric_SSIM, 'LineWidth', 2);
    title('SSIM (maximize)', 'FontSize', 18);
    xticks(min(accRange):gap:max(accRange));
    set(gca, 'LineWidth', 1.5);
    
    % % Plot JS Divergence.
    % subplot(4, 1, 4)
    % plot(accRange, Metric_JSD, 'LineWidth', 2);
    % title('Log(JS Divergence) (minimize)', 'FontSize', 18);
    % xticks(min(accRange):gap:max(accRange));
    % set(gca, 'LineWidth', 1.5);
end

%% ----- Combine Metrics into Output Matrix -----
% Each row corresponds to an element in accRange with columns:
% [L2 distance, PCC, Cosine similarity]
Metric_optz = cat(2, Metric_distance, Metric_pcc, Metric_SSIM);

>>>>>>> 5e4cf7b (Initial commit)
end