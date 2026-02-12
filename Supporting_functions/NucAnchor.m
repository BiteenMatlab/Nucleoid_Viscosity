function [anchor_modi, kb] = NucAnchor(outline, aspt_ratio, convert_factor, num_pix_y, fig_plot)
% NucAnchor Processes a nucleoid outline to extract “anchor dots” by
% removing redundant collinear points. The resulting anchor dots are symmetrized
% across the vertical axis.
%
% [anchor_modi, kb] = NucAnchor(outline, aspt_ratio, convert_factor, num_pix_y, fig_plot)
%
% Inputs:
% outline - An [N x 2] numeric array of (x,y) coordinates of the nucleoid outline.
% aspt_ratio - Scalar representing the aspect ratio (e.g., 2).
% convert_factor- Conversion factor to scale the outline coordinates.
% num_pix_y - (Optional) Scalar (default 20) representing the number of pixels in y–direction.
% fig_plot - (Optional) String controlling plotting. Options:
% 'original' - Plot raw (shifted/scaled) outline dots.
% 'anchor' - Plot anchor dots computed by AnchorDot.
% 'anchor_final' - Plot the final modified anchor dots.
% 'off' - No plotting.
%
% Outputs:
% anchor_modi - An M x 2 array of processed (anchor) coordinates (symmetrized).
% kb - An (M-1) x 2 array of line properties (slope and intercept) between consecutive anchor dots.
%
% Example:
% [anchor_modi, kb] = NucAnchor(nuc_outline, 2, 1, 20, 'anchor_final');
%
% Author: Xiaofeng Dai
% Date: 04/15/2025

%% 1. Set Defaults for Optional Inputs
if nargin == 3
    num_pix_y = 20;
    fig_plot = 'off';
elseif nargin == 4
    fig_plot = 'off';
end

%% 2. Shift and Scale the Outline
% Calculate shift values to center the outline.
shiftx = num_pix_y * aspt_ratio / 2 + 0.5;
shifty = num_pix_y / 2 + 0.5;
outline(:,1) = outline(:,1) - shiftx;
outline(:,2) = outline(:,2) - shifty;

% Keep only points on or above the horizontal axis.
outline(outline(:,2) < 0, :) = [];

% Flip the outline such that the order is reversed.
outline = flip(outline);

% Scale the outline using the conversion factor.
outline = outline * convert_factor;

% Ensure that the outline intersects the horizontal axis.
buff_ind = find(outline(:,2)==0 & outline(:,1) < 0, 1);
if isempty(buff_ind)
    % If no intersection is found, create two points at the start and end.
    buff1 = [outline(1,1), 0];
    buff2 = [outline(end,1), 0];
    outline = cat(1, buff1, outline, buff2);
end
clear buff_ind buff1 buff2 shiftx shifty

%% 3. Optional Plot: Show the "Original" Outline
if strcmpi(fig_plot, 'original')
    figure;
    scatter(outline(:,1), outline(:,2), 5, 'filled');
    axis equal;
    xlim([min(outline(:,1))-2, max(outline(:,1))+2]);
    ylim([min(outline(:,2))-1, max(outline(:,2))+1]);
    title([fig_plot, ' dots from nucleoid outline']);
end

%% 4. Compute Anchor Dots using AnchorDot
[anchor_coord, slope] = AnchorDot(outline);

if strcmpi(fig_plot, 'anchor')
    figure;
    scatter(outline(:,1), outline(:,2), 5, 'filled');
    hold on;
    scatter(anchor_coord(:,1), anchor_coord(:,2), 15, 'red', 'x');
    plot(anchor_coord(:,1), anchor_coord(:,2), 'r');
    axis equal;
    xlim([min(outline(:,1))-2, max(outline(:,1))+2]);
    ylim([min(outline(:,2))-1, max(outline(:,2))+1]);
    title([fig_plot, ' dots']);
end

%% 5. Adjust Anchor Coordinates for Vertical Segments
% For segments with near-vertical slopes (abs(slope) > 1e8), nudge the x coordinate slightly
% and recompute slopes until no such segments remain.
while ~isempty(find(abs(slope) > 1e8, 1))
    for i = 1:length(slope)
        if abs(slope(i)) > 1e8  % indicates a vertical line segment
            anchor_coord(i,1) = anchor_coord(i,1) - 1e-10;
        end
    end
    % Recalculate slopes for the modified anchors.
    tt = diff(anchor_coord);
    tt(tt(:,1)==0,1) = 1e-15;
    slope = tt(:,2) ./ tt(:,1);
end
clear tt i

%% 6. Generate Symmetric Anchor Outline
% Extract anchors with negative x-coordinates.
left_anchors = anchor_coord(anchor_coord(:,1) < 0, :);
% Sort the left anchors in ascending order (by x, then y).
left_anchors = sortrows(left_anchors);
% Mirror the left anchors by changing sign of x-coordinate.
right_anchors = left_anchors;
right_anchors(:,1) = -right_anchors(:,1);
right_anchors = flip(right_anchors);
% Concatenate left and mirrored right anchors.
anchor_modi = [left_anchors; right_anchors];
clear left_anchors right_anchors

% remove duplicate points
[C, ~, ic] = unique(anchor_modi(:,1)); 
counts = histcounts(ic, 1:max(ic)+1);
duplicates = C(counts > 1);
for i = 1:length(duplicates)
    dup = duplicates(i);
    dup_ind = find(anchor_modi(:,1)==dup);
    buff = mean(anchor_modi(dup_ind,2));
    anchor_modi(dup_ind,1) = nan;
    anchor_modi(dup_ind(1),1) = dup;
    anchor_modi(dup_ind(1),2) = buff;
end
aa = isnan(anchor_modi(:,1));
anchor_modi(aa,:) = [];

%%% sort anchor_modi by minimal adjacent distance
% comment out if not stable
new_order = anchor_modi(find(anchor_modi(:,2) == 0, 1, 'first'),:);
dist_store = zeros(height(anchor_modi)-1,1);
copy_anchor_modi = anchor_modi;
for iter_dis = 1:(height(anchor_modi)-1)
    now_dot = new_order(end,:);
    del_row = ismember(copy_anchor_modi,now_dot, 'rows');
    copy_anchor_modi(del_row,:) = [];
    dis_list = sum((now_dot - copy_anchor_modi).^2,2);
    next_rowid = dis_list == min(dis_list);
    dist_store(iter_dis) = dis_list(next_rowid);
    next_dot = copy_anchor_modi(next_rowid,:);
    new_order = cat(1,new_order,next_dot);
end
new_order(find(new_order(:,2) == 0, 1, 'last')+1:end,:)=[];
anchor_modi = new_order;
%%%

%% 7. Optional Plot: Display Final Anchor Dots
if strcmpi(fig_plot, 'anchor_final')
    figure;
    scatter(outline(:,1), outline(:,2), 5, 'filled');
    hold on;
    scatter(anchor_modi(:,1), anchor_modi(:,2), 15, 'red', 'x');
    plot(anchor_modi(:,1), anchor_modi(:,2), 'r-', 'LineWidth', 1.2);
    axis equal;
    xlim([min(outline(:,1))-2, max(outline(:,1))+2]);
    ylim([min(outline(:,2))-1, max(outline(:,2))+1]);
    title('Final Anchor dots');
end

%% 8. Compute Line Properties from Anchor Dots
kb = LineProp(anchor_modi);
% Remove any segments resulting in NaN values.
nanIdx = find(isnan(kb));
if ~isempty(nanIdx)
    anchor_modi(nanIdx, :) = [];
    kb = LineProp(anchor_modi);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kb = LineProp(anchor)
% LineProp Computes the slope and intercept for each segment connecting consecutive anchor points.
%
% kb = LineProp(anchor) returns an (M-1) x 2 array where each row contains:
% [slope, intercept]
% computed from consecutive points in the anchor array.
%
% Input:
% anchor - An M x 2 numeric array of anchor coordinates.
%
% Output:
% kb - An (M-1) x 2 numeric array such that:
% kb(i,1): slope between anchor(i,:) and anchor(i+1,:)
% kb(i,2): intercept of the line segment.
%
% Author: Xiaofeng Dai
% Date: 04/15/2025

N = size(anchor, 1);
kb = zeros(N-1, 2);
for i = 1:(N - 1)
    pt1 = anchor(i, :);
    pt2 = anchor(i + 1, :);
    % Avoid division by zero by substituting a small number if necessary.
    dx = pt2(1) - pt1(1);
    if abs(dx) < 1e-10
        dx = 1e-10;
    end
    kb(i, 1) = (pt2(2) - pt1(2)) / dx;
    % Compute intercept using: b = y - slope*x.
    kb(i, 2) = pt1(2) - kb(i, 1) * pt1(1);
end
end