function [avg_aspect_ratio, avg_width, cellIndices] = cellMorphFilter(sel_filter, bounds, cell_morph, morph_by_mov, num_pix_y)
% cellMorphFilter Filters cells based on a specified morphology parameter and
% computes the average aspect ratio for cells that pass the filter.
%
% [avg_aspect_ratio, cellIndices] = cellMorphFilter(sel_filter, bounds, cell_morph, morph_by_mov)
%
% INPUTS:
% sel_filter - Scalar specifying which column of cell_morph (and morph_by_mov)
% is used for filtering. For example:
% 2 --> cell length,
% 3 --> cell width,
% 4 --> aspect ratio,
% 5 --> cell area.
% bounds - 1x2 numeric vector [low_bd, up_bd]. Only cells with values in
% cell_morph(:,sel_filter) strictly between low_bd and up_bd pass.
% cell_morph - Numeric matrix containing morphology data for all cells.
% It is assumed that column 4 holds the aspect ratio.
% morph_by_mov- Cell array. Each element is a matrix with morphology metrics
% (with the first column being the cell ID) for cells from one ROI.
% num_pix_y - (Optional) Scalar specifying number of grids in y–axis.
% Default is 20.
%
% OUTPUTS:
% avg_aspect_ratio - The average aspect ratio (from column 4 of cell_morph)
% avg_width - The average width (from column 3 of cell_morph)
% for cells passing the filter. It is quantized based on num_pix_y.
% cellIndices - A cell array (one element per ROI) containing the indices
% (cell IDs from column 1) of cells that pass the filter.
%
% The function also prints to the command window the mean aspect ratio and the number
% of grids in both x- and y–axes.
%
% Example:
% % Filter cells based on cell length between 1.2 and 2.2:
% [avg_ar, idx] = cellMorphFilter(2, [1.2, 2.2], cell_morph, morph_by_mov);
%
% Author: Xiaofeng Dai
% Date: 04/09/2025

%% 1. Input Validation and Default
if nargin < 4
    error('At least four inputs are required: sel_filter, bounds, cell_morph, and morph_by_mov');
end
if nargin < 5 || isempty(num_pix_y)
    num_pix_y = 20;
end
if ~isscalar(sel_filter) || ~ismember(sel_filter, [2,3,4,5])
    error('sel_filter must be one of the following scalar values: 2, 3, 4, or 5.');
end
if numel(bounds) ~= 2
    error('bounds must be a 1x2 vector, e.g. [low_bd, up_bd].');
end
low_bd = bounds(1);
up_bd  = bounds(2);

if ~isnumeric(cell_morph) || isempty(cell_morph)
    error('cell_morph must be a non-empty numeric matrix.');
end
if ~iscell(morph_by_mov) || isempty(morph_by_mov)
    error('morph_by_mov must be a non-empty cell array.');
end

%% 2. Compute Averaged Aspect Ratio and Width for Cells Passing the Filter
% Use the chosen column from cell_morph to determine which cells pass.
filterValues = cell_morph(:, sel_filter);
passFilter   = filterValues > low_bd & filterValues < up_bd;

% Aspect ratio is stored in column 4.
mean_ar = mean(cell_morph(passFilter, 4));
% Width is stored in column 3.
avg_width = mean(cell_morph(passFilter, 3));

% Quantize the mean aspect ratio to an integer number of grids along y axis.
avg_aspect_ratio = round(mean_ar * num_pix_y) / num_pix_y;

% Display the results.
fprintf('Mean aspect ratio = %.2f    ~%.2f\n', mean_ar, avg_aspect_ratio);
fprintf('Number of grids y axis = %d\nNumber of grids x axis = %.0f\n', num_pix_y, avg_aspect_ratio*num_pix_y);

%% 3. Determine the Indices of Cells That Pass the Filter for Each ROI
% The number of movies is determined from the length of morph_by_mov.
num_movies = numel(morph_by_mov);
cellIndices = cell(num_movies, 1);

for i = 1:num_movies
    roiMorph = morph_by_mov{i};
    % Determine the cells in this movie that pass the filter.
    % test1: values less than up_bd, test2: greater than low_bd.
    test1 = roiMorph(:, sel_filter) < up_bd;
    test2 = roiMorph(:, sel_filter) > low_bd;
    passROI = (test1 & test2);  % Logical AND instead of adding tests.
    
    % Assume that the cell ID is stored in the first column.
    roiIndices = roiMorph(passROI, 1);
    cellIndices{i} = roiIndices;
end
end