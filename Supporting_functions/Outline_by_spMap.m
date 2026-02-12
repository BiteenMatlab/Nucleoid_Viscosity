function [nuc_contour, up_bd, low_bd, g] = Outline_by_spMap(spMap, cntr_range, figplot, cell_outline, aspt_ratio, num_pix_y)
% Outline_by_spMap Calculates nucleoid boundaries using a spatial displacement map.
%
% [nuc_contour, up_bd, low_bd, g] = Outline_by_spMap(spMap, cntr_range, figplot,
% cell_outline, aspt_ratio, num_pix_y) performs an interpolation and smoothing
% of the spatial displacement map (spMap) and, using either a provided 1×2 vector
% (cntr_range) or coordinates (from htMap analysis), sets upper (up_bd) and lower
% (low_bd) bounds for contour thresholding. Contour lines are then computed from the
% gradient of the smoothed map; these outlines are returned in nuc_contour. 
% If figplot is 'on', the original spMap is displayed with the nucleoid boundaries overlaid.
%
% INPUTS:
% spMap - [m×n] numeric matrix representing the original spatial displacement map.
% cntr_range - Either a 1×2 vector with direct upper and lower bounds (in normalized units)
% or a set of coordinates (e.g. from htMap analysis) used to calculate the bounds.
% figplot - (Optional) 'on' (default) or 'off'; if 'on' the function plots the outline.
% cell_outline - (Optional) [P×2] numeric matrix with cell contour coordinates to overlay.
% aspt_ratio - (Optional) Scalar (default 2) for x–axis scaling relative to y.
% num_pix_y - (Optional) Scalar (default 20) giving the number of grid bins in the y–direction.
%
% OUTPUTS:
% nuc_contour - A 2×1 cell array. nuc_contour{1} and nuc_contour{2} contain the (x,y)
% coordinates (normalized) of two clusters of contour points.
% up_bd - Upper bound (normalized, between 0 and 1) used for contour thresholding.
% low_bd - Lower bound used for contour thresholding.
% g - Figure handle if a plot was produced (empty otherwise).
%
% Example:
% [nuc_outline, up_bd, low_bd, g] = Outline_by_spMap(spMap, [0.9 0.8], 'on', cell_outline, 2, 20);
%
% Author: Xiaofeng Dai
% Date: 04/12/2025

%% 1. Set Defaults and Validate Inputs
if nargin < 2
    error('At least two inputs (spMap and cntr_range) are required.');
end
if nargin < 3 || isempty(figplot),    figplot = 'on';    end
if nargin < 4,                         cell_outline = []; end
if nargin < 5 || isempty(aspt_ratio),   aspt_ratio = 2;    end
if nargin < 6 || isempty(num_pix_y),    num_pix_y = 20;    end
if ~isnumeric(spMap)
    error('spMap must be a numeric matrix.');
end

%% 2. Interpolate & Smooth the Spatial Displacement Map
% Replace NaNs with zeros.
spMap_clean = spMap;
spMap_clean(isnan(spMap_clean)) = 0;
% Interpolate spMap to a grid that is 10× finer using cubic interpolation.
sp_intp = Map_Interpolate(spMap_clean, 10, 'cubic');
% Smooth the interpolated map using a Gaussian filter.
sp_intp_smth = smoothdata2(sp_intp, "gaussian", 'SmoothingFactor', 0.1);
% Scale the smoothed map back (since Map_Interpolate does not change total magnitude):
sp_intp_smth = sp_intp_smth * 100;  % equivalent to 10*10

%% 3. Determine Threshold Bounds from cntr_range
% Two cases: if cntr_range is a 1×2 vector, use its max/min directly;
% otherwise, assume it contains coordinate data from which bounds are computed.
if ~isequal(size(cntr_range), [1,2])
    % Convert normalized coordinates to index space in sp_intp_smth.
    a = (cntr_range - 0.5) * 10;  % scale factor matching the interpolation above
    numPts = size(a,1);
    b = zeros(numPts,1);
    for i = 1:numPts
        % Round indices to ensure valid indexing.
        rowIdx = max(1, min(round(a(i,2)), size(sp_intp_smth,1)));
        colIdx = max(1, min(round(a(i,1)), size(sp_intp_smth,2)));
        b(i) = sp_intp_smth(rowIdx, colIdx);
    end
    % Normalize b by the maximum value of the smoothed map.
    b = b / max(sp_intp_smth(:));
    % Set the threshold bounds (ensure within [0,1]).
    up_bd = min(ceil(max(b)*100)/100, 1);
    low_bd = max(floor(min(b)*100)/100, 0);
else
    up_bd = max(cntr_range);
    low_bd = min(cntr_range);
end

%% 4. Compute Contour from the Gradient of the Smoothed Map
% Set contour levels based on the thresholds scaled by the maximum intensity.
cnt_level2 = [up_bd, low_bd] * max(sp_intp_smth(:));
% Compute contour lines.
M1 = contourc(sp_intp_smth, cnt_level2);
% Rescale contour coordinates to normalized [0,1] given the interpolation procedure.
M1 = M1 ./ 10 + 0.5;
M1 = M1.';  % transpose so each row is an (x,y) pair.
% Cluster the contour points using DBSCAN (epsilon = 0.5, minimum points = 2).
clusterIdx = dbscan(M1, 0.5, 2);
M11 = M1(clusterIdx == 1, :);
M12 = M1(clusterIdx == 2, :);
nuc_contour = cell(2, 1);
nuc_contour{1} = M11;
nuc_contour{2} = M12;

%% 5. Plot (if Requested)
if strcmpi(figplot, 'on')
    g = figure;
    imagesc(spMap);    % Plot original spatial displacement map.
    view([0 90]);
    hold on;
    % Overlay the two nucleoid contour clusters with different line styles.
    plot(M11(:,1), M11(:,2), 'k-', 'LineWidth', 1.2);
    plot(M12(:,1), M12(:,2), 'k--', 'LineWidth', 1.2);
    view([0 90]);
    pbaspect([aspt_ratio, 1, 1]);
    hold on;
    % If a cell outline was provided, overlay it.
    if ~isempty(cell_outline)
        plot(cell_outline(:,1), cell_outline(:,2), 'white', 'LineWidth', 1.5);
        % Connect first and last point.
        plot([cell_outline(end,1); cell_outline(1,1)], [cell_outline(end,2); cell_outline(1,2)], 'white', 'LineWidth', 1.5);
    end
    % Set colormap and color limits.
    colormap(slanCM('parula'));
    clim([10 * floor(min(spMap(:))/10), 10 * ceil(max(spMap(:))/10)]);
    % Define plot limits and tick labels.
    num_pix_x = num_pix_y * aspt_ratio;
    xlim([0.5, num_pix_x + 0.5]);
    xticks(linspace(0.5, num_pix_x + 0.5, 6));
    ylim([0.5, num_pix_y + 0.5]);
    yticks(linspace(0.5, num_pix_y + 0.5, 6));
    xticklabels({'0','0.2','0.4','0.6','0.8','1'});
    yticklabels({'1','0.8','0.6','0.4','0.2','0'});
    colorbar;
else
    g = [];
end
end