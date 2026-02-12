function fig_htMap(mapData, numPixelsY, aspectRatio, colorMap, showColorbar, colorRange)
% fig_htMap displays a heat map with axes formatting and optional color limits.
%
% Inputs:
%   mapData      - Numeric matrix for the heat map.
%   numPixelsY   - Number of pixels along the y-axis.
%   aspectRatio  - Aspect ratio for displaying the heat map.
%   colorMap     - Colormap for the display.
%   showColorbar - 'on' to display a colorbar; 'off' otherwise.
%   colorRange   - Two-element vector [min max] for color limits. If empty, no limits are applied.
%
% This function uses imagesc to display the heat map, sets axis ticks/labels,
% applies the aspect ratio, and optionally imposes the color limits.
% Author: Xiaofeng Dai
% Date: 05/18/2025

% Determine the number of pixels along the x-axis.
numPixelsX = numPixelsY * aspectRatio;

% Display the map.
imagesc(mapData);
axis equal;
colormap(gca, colorMap);

% Display the colorbar if requested.
if strcmpi(showColorbar, 'on')
    colorbar;
end

% Set x-axis limits and ticks.
xlim([0.5, numPixelsX + 0.5])
xticks(linspace(0.5, numPixelsX + 0.5, 6))

% Set y-axis limits and ticks.
ylim([0.5, numPixelsY + 0.5])
yticks(linspace(0.5, numPixelsY + 0.5, 6))

% Define custom tick labels.
xticklabels({'0','0.2','0.4','0.6','0.8','1'});
yticklabels({'1','0.8','0.6','0.4','0.2','0'});

% Adjust the plot aspect ratio.
pbaspect([aspectRatio 1 1])

% Apply color limits if provided.
if ~isempty(colorRange)
    clim(gca, colorRange);
end
end