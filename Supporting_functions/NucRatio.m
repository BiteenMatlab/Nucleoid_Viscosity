function [nuc_ratio, cell_para, nuc_para] = NucRatio(cell_bndy, nuc_bndy, aspt_ratio, num_pix_y, cell_rodfit)
% NucRatio Computes relative morphology ratios between nucleoid and cell.
%
% [nuc_ratio, cell_para, nuc_para] = NucRatio(cell_bndy, nuc_bndy, aspt_ratio, num_pix_y, cell_rodfit)
%
% This function calculates the relative length and width ratios (as well as the area
% ratio) between the nucleoid and cell boundary. Two methods for determining cell
% dimensions are available via the 'cell_rodfit' flag:
%
% 'on' - Uses a rod‐model (via Rod_Morph_fit) to estimate cell boundaries. The cell
% length is computed as 2*(R + L) and the width as 2*R.
%
% 'off' - Uses direct measurements from the boundary coordinates.
%
% Next, the function computes the nucleoid length and width from nuc_bndy. Both cell
% parameters and nucleoid parameters are then normalized using the cell width. Finally,
% the ratios (length and width) are computed from the normalized values, and an area
% ratio (nucleoid area divided by cell area) is appended.
%
% Inputs:
% cell_bndy - [N x 2] numeric array of cell boundary coordinates.
% nuc_bndy - [M x 2] numeric array of nucleoid boundary coordinates.
% aspt_ratio - Scalar specifying the cell aspect ratio (x/y).
% num_pix_y - (Optional) Scalar, e.g., 20 (used for normalization).
% cell_rodfit- (Optional) 'on' (default) to perform rod‐fitting; 'off' to use direct measurements.
%
% Outputs:
% nuc_ratio - 1×3 vector: [nuc_length_ratio, nuc_width_ratio, nuc_area_ratio]
% cell_para - 1×2 vector: [cell_length, cell_width] normalized by cell width.
% nuc_para - 1×2 vector: [nuc_length, nuc_width] normalized by cell width.
%
% Example:
% [nuc_ratio, cell_para, nuc_para] = NucRatio(cell_bndy, nuc_bndy, 2);
%
% Author: Xiaofeng Dai
% Date: 04/15/2025

%% 1. Set Default Values for Optional Inputs
if nargin < 4 || isempty(num_pix_y)
    num_pix_y = 20;
end
if nargin < 5 || isempty(cell_rodfit)
    cell_rodfit = 'on';
end

%% 2. Compute Cell Parameters
switch lower(cell_rodfit)
    case 'on'
        % Use rod fitting method.
        [Rfit, Lfit] = Rod_Morph_fit(cell_bndy, num_pix_y, aspt_ratio, 'off');
        cell_length = 2 * (Rfit + Lfit); % Total cell length
        cell_width  = 2 * Rfit;           % Total cell width
    case 'off'
        % Use direct measurement: cell length from x-span.
        cell_length = max(cell_bndy(:,1)) - min(cell_bndy(:,1));
        % Compute cell width as the difference between the mean of top 10% and bottom 10% y-values.
        pct = round(size(cell_bndy, 1) * 0.1);
        cell_width = mean(maxk(cell_bndy(:,2), pct)) - mean(mink(cell_bndy(:,2), pct));
    otherwise
        error('cell_rodfit must be either ''on'' or ''off''.');
end
cell_para = [cell_length, cell_width];

%% 3. Compute Nucleoid Parameters from nucleoid boundary
nuc_length = max(nuc_bndy(:,1)) - min(nuc_bndy(:,1));
pct = round(size(nuc_bndy, 1) * 0.1);
nuc_width = mean(maxk(nuc_bndy(:,2), pct)) - mean(mink(nuc_bndy(:,2), pct));
nuc_para = [nuc_length, nuc_width];

%% 4. Normalize Parameters by Cell Width
cell_para = cell_para ./ cell_width;  % Normalize cell dimensions
nuc_para  = nuc_para ./ cell_width;    % Normalize nucleoid dimensions

%% 5. Compute Ratio Metrics
% Relative length and width ratios.
nuc_ratio = nuc_para ./ cell_para;
% Compute relative area ratio.
area_cell = polyarea(cell_bndy(:,1), cell_bndy(:,2));
area_nuc  = polyarea(nuc_bndy(:,1), nuc_bndy(:,2));
area_ratio = area_nuc / area_cell;
% Append area ratio.
nuc_ratio = [nuc_ratio, area_ratio];
end
