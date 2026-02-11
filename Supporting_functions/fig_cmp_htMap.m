function fig_cmp_htMap(expHeatMap, simuHeatMaps, accessibilityRange, selectedAcc, numPixelsY, aspectRatio, colorMap, singleColorBar, colorRange)
% fig_cmp_htMap displays experimental and simulated heat maps side-by-side.
%
%   fig_cmp_htMap(expHeatMap, simuHeatMaps, accessibilityRange, selectedAcc, 
%                 numPixelsY, aspectRatio, colorMap, singleColorBar, colorRange)
%
% Inputs:
%   expHeatMap         - Matrix representing the experimental heat map data.
%   simuHeatMaps       - Cell array of simulated heat map matrices. Each cell 
%                        corresponds to an accessibility value.
%   accessibilityRange - Numeric vector of accessibility values for simuHeatMaps.
%   selectedAcc        - Numeric vector of selected accessibility values that
%                        you want to display.
%   numPixelsY         - Number of pixels along the Y-axis (for heat map display).
%   aspectRatio        - Aspect ratio for the heat map.
%   colorMap           - Colormap used for displaying the heat maps.
%   singleColorBar     - Logical flag; if true, a common colorbar is used for all subplots.
%   colorRange         - Two-element vector [min max] specifying color limits.
%                        If empty and singleColorBar is true, the limits are determined from data.
%
% This function generates a figure with subplots. The first subplot always displays
% the experimental heat map, and the remaining subplots display simulated heat maps which
% correspond to the selected accessibility values.
%
% Author: Xiaofeng Dai
% Date: 05/18/2025

%% ----- Defensive Checks -----
% Check that expHeatMap is numeric.
if ~isnumeric(expHeatMap)
    error('expHeatMap must be a numeric matrix.');
end
% Check that simuHeatMaps is a cell array.
if ~iscell(simuHeatMaps)
    error('simuHeatMaps must be a cell array.');
end
% Check that accessibilityRange is numeric and that it matches the size of simuHeatMaps.
if ~isnumeric(accessibilityRange) || numel(accessibilityRange) ~= numel(simuHeatMaps)
    error('accessibilityRange must be a numeric vector with the same number of elements as simuHeatMaps.');
end
% Check that selectedAcc is numeric.
if ~isnumeric(selectedAcc)
    error('selectedAcc must be a numeric vector.');
end
% Check numPixelsY and aspectRatio.
if ~isscalar(numPixelsY) || numPixelsY <= 0
    error('numPixelsY must be a positive scalar.');
end
if ~isscalar(aspectRatio) || aspectRatio <= 0
    error('aspectRatio must be a positive scalar.');
end
% Check singleColorBar is a logical value.
if ~islogical(singleColorBar)
    error('singleColorBar must be a logical value (true or false).');
end
% If provided, colorRange must be a two-element numeric vector.
if ~isempty(colorRange) && (~isnumeric(colorRange) || numel(colorRange) ~= 2)
    error('colorRange must be a 2-element numeric vector if provided.');
end

%% ----- Prepare Simulated Heat Maps for Selected Accessibility Values -----
% Allocate a cell array to store the simulated heat maps for the selected accessibility values.
numSelectedAcc = numel(selectedAcc);
cmpAccHeatMaps = cell(numSelectedAcc,1);

for i = 1:numSelectedAcc
    % Current selected accessibility value.
    accValue = selectedAcc(i);
    % Identify the corresponding index from accessibilityRange (using a tolerance).
    idxLogical = abs(accessibilityRange - accValue) < 1e-10;
    if sum(idxLogical) ~= 1
        error('Each selected accessibility value must match exactly one element in accessibilityRange.');
    end
    % Extract the simulated heat map corresponding to this accessibility value.
    cmpAccHeatMaps{i} = simuHeatMaps{idxLogical};
end

%% ----- Determine Color Limits (if using a single colorbar) -----
if singleColorBar
    if isempty(colorRange)
        % Combine data from each simulated heat map and the experimental heat map.
        combinedData = cat(1, cmpAccHeatMaps{:}, expHeatMap);
        minColorLimit = min(combinedData, [], 'all');
        maxColorLimit = max(combinedData, [], 'all');
    else
        minColorLimit = colorRange(1);
        maxColorLimit = colorRange(2);
    end
end

%% ----- Plotting: Experimental and Simulated Heat Maps -----
f = figure;
numSubplots = numSelectedAcc + 1;   % one subplot for experimental data and one per selected simulation
numRows = ceil(numSubplots / 2);      % Arrange subplots in two columns

% Plot the Experimental Heat Map in the first subplot.
subplot(numRows, 2, 1)
if singleColorBar
    fig_htMap(expHeatMap, numPixelsY, aspectRatio, colorMap, 'off', [minColorLimit, maxColorLimit])
else
    fig_htMap(expHeatMap, numPixelsY, aspectRatio, colorMap, 'on', [])
end
title('Experimental Data', 'FontWeight', 'bold','FontSize',14)
hold on

% Plot each selected Simulated Heat Map.
for i = 1:numSelectedAcc
    subplot(numRows, 2, i+1)
    if singleColorBar
        fig_htMap(cmpAccHeatMaps{i}, numPixelsY, aspectRatio, colorMap, 'off', [minColorLimit, maxColorLimit])
    else
        fig_htMap(cmpAccHeatMaps{i}, numPixelsY, aspectRatio, colorMap, 'on', [])
    end
    title({['Accessibility = ' sprintf('%.2f', selectedAcc(i))]}, 'FontWeight', 'bold','FontSize',14)
    hold on
end

% Set the overall figure position.
f.Position = [100 100 800 numRows*200+50];

% If a single colorbar is used, attach it to the current axes.
if singleColorBar
    ax = gca;
    cBar = colorbar(ax, 'Position', [0.92 0.1 0.03 0.82]);
    colormap(cBar, colorMap)
    clim(ax, [minColorLimit, maxColorLimit]);
end
end