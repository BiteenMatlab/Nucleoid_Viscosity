function htMap_final_bndy(Map, cell_outline, nuc_contour, num_pix_y, aspt_ratio, clr_map)
% htMap_final_bndy Plots a heat map with overlaid nucleoid and cell outlines.
%
% htMap_final_bndy(Map, cell_outline, nuc_contour, num_pix_y, aspt_ratio, clr_map)
%
% This function displays the input heat map (Map) and overlays two outlines:
% - nuc_contour: the nucleoid outline drawn in red dashed lines.
% - cell_outline: the cell outline drawn in white solid lines.
% The axes are scaled based on the given number of y-grid bins (num_pix_y) and
% aspect ratio (aspt_ratio). A specified colormap (clr_map) is applied.
%
% INPUTS:
% Map - A numeric matrix representing the heat map.
% cell_outline- An [N x 2] numeric array containing (x,y) cell outline coordinates.
% nuc_contour- An [M x 2] numeric array containing (x,y) nucleoid contour coordinates.
% num_pix_y - Scalar; number of bins along the y-axis (e.g., 20).
% aspt_ratio - Scalar; ratio for the x-axis relative to y (e.g., 2).
% clr_map - A valid MATLAB colormap (e.g., 'jet' or a custom colormap).
%
% Example:
% htMap_final_bndy(heatMapMatrix, cellOutline, nucContour, 20, 2, 'jet');
%
% Author: Xiaofeng Dai
% Date: 04/13/2025

%% 1. Input Validation
if ~ismatrix(Map) || ~isnumeric(Map)
    error('Map must be a numeric 2D matrix.');
end
if size(cell_outline,2) ~= 2 || size(nuc_contour,2) ~= 2
    error('cell_outline and nuc_contour must be [N x 2] numeric arrays.');
end
if ~isscalar(num_pix_y) || ~isscalar(aspt_ratio)
    error('num_pix_y and aspt_ratio must be scalars.');
end

%% 2. Create Figure and Display Heat Map
figure;
imagesc(Map);
pbaspect([aspt_ratio, 1, 1]);
hold on;

%% 3. Overlay Outlines
% Plot nucleoid outline (red dashed lines).
plot(nuc_contour(:,1), nuc_contour(:,2), 'r--', 'LineWidth', 1.2);
% Plot cell outline (white solid lines); also connect the last point to the first.
plot(cell_outline(:,1), cell_outline(:,2), 'w-', 'LineWidth', 1.2);
plot([cell_outline(end,1); cell_outline(1,1)], [cell_outline(end,2); cell_outline(1,2)],...
     'w-', 'LineWidth', 1.2);

%% 4. Apply Colormap and Set Axis Properties
colormap(clr_map);
num_pix_x = num_pix_y * aspt_ratio;
xlim([0.5, num_pix_x + 0.5]);
xticks(linspace(0.5, num_pix_x + 0.5, 6));
ylim([0.5, num_pix_y + 0.5]);
yticks(linspace(0.5, num_pix_y + 0.5, 6));
xticklabels({'0','0.2','0.4','0.6','0.8','1'});
yticklabels({'1','0.8','0.6','0.4','0.2','0'});
colorbar;

hold off;
end