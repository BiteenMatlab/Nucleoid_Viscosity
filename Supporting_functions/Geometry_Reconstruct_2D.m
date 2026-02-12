function Geometry_Reconstruct_2D(Radius, Length, scalingCellWidth, scalingCellLength, anchor_htMap, anchor_spMap)
% GEOMETRY_RECONSTRUCT_2D Plots 2D representations of the cell boundary,
% inner membrane, and nucleoid boundaries.
%
%   Geometry_Reconstruct_2D(Radius, Length, scalingCellWidth, scalingCellLength, anchor_htMap, anchor_spMap)
%
%   Input arguments:
%     Radius            - Effective radius for the inner membrane boundary.
%     Length            - Half-length of the cylinder part of inner membrane boundary.
%     scalingCellWidth  - Parameter defining the cell’s width (used for the outer boundary).
%     scalingCellLength - Parameter defining the cell’s overall length (used for the outer boundary).
%     anchor_htMap      - A [n x 2] anchor map for the nucleoid boundary from the heat map.
%                         
%     anchor_spMap      - A [m x 2] anchor map for the nucleoid boundary from the spatial displacement map.
%                         Same format as anchor_htMap.
%
%
% Author: Xiaofeng Dai
% Date: 05/10/2025

    %% --- Defensive input validation ---
    narginchk(6,6);

    % Validate anchor maps: check that they have at least two columns and two rows.
    if size(anchor_htMap,2) < 2 || size(anchor_spMap,2) < 2
        error('Both anchor_htMap and anchor_spMap must have at least two columns [independent, dependent].');
    end
    if size(anchor_htMap,1) < 2 || size(anchor_spMap,1) < 2
        error('Both anchor_htMap and anchor_spMap must have at least two rows.');
    end

    % Ensure that the anchor maps are sorted by their independent variable (first column).
    if ~issorted(anchor_htMap(:,1))
        warning('anchor_htMap first column is not sorted. Sorting anchor_htMap now.');
        anchor_htMap = sortrows(anchor_htMap,1);
    end
    if ~issorted(anchor_spMap(:,1))
        warning('anchor_spMap first column is not sorted. Sorting anchor_spMap now.');
        anchor_spMap = sortrows(anchor_spMap,1);
    end

    %% --- Define boundary functions ---

    % Outer cell boundary (using scaling parameters)
    % The outer cell boundary is defined by a custom function "set_cellx".
    outerCellBoundary = @(y) set_cellx(scalingCellWidth, scalingCellLength - scalingCellWidth, y);
    y_outer = linspace(-scalingCellLength, scalingCellLength, 100);
    x_outer = outerCellBoundary(y_outer);

    % Inner membrane boundary (using Radius and Length parameters)
    innerBoundary = @(y) set_cellx(Radius, Length, y);
    y_inner = linspace(-(Radius + Length), (Radius + Length), 100);
    x_inner = innerBoundary(y_inner);

    % Nucleoid boundary from heat map
    nucBoundary_ht = @(y) piecefunc(anchor_htMap, y);
    y_nuc_ht = linspace(min(anchor_htMap(:,1)), max(anchor_htMap(:,1)), 100);
    x_nuc_ht = nucBoundary_ht(y_nuc_ht);

    % Nucleoid boundary from spatial displacement map
    nucBoundary_sp = @(y) piecefunc(anchor_spMap, y);
    y_nuc_sp = linspace(min(anchor_spMap(:,1)), max(anchor_spMap(:,1)), 100);
    x_nuc_sp = nucBoundary_sp(y_nuc_sp);

    %% --- Visualization ---
    figure;
    hold on;

    % Plot outer cell boundary (both positive and negative x)
    plot(y_outer, x_outer, 'Color', "#0072BD", 'LineWidth', 1, 'LineStyle', '--');
    plot(y_outer, -x_outer, 'Color', "#0072BD", 'LineWidth', 1, 'LineStyle', '--');

    % Plot inner membrane boundary
    plot(y_inner, x_inner, 'Color', "#4DBEEE", 'LineWidth', 1.5, 'LineStyle', '-');
    plot(y_inner, -x_inner, 'Color', "#4DBEEE", 'LineWidth', 1.5, 'LineStyle', '-');

    % Plot nucleoid boundary from height map
    plot(y_nuc_ht, x_nuc_ht, 'Color', "#A2142F", 'LineWidth', 1, 'LineStyle', ':');
    plot(y_nuc_ht, -x_nuc_ht, 'Color', "#A2142F", 'LineWidth', 1, 'LineStyle', ':');

    % Plot nucleoid boundary from spot map
    plot(y_nuc_sp, x_nuc_sp, 'Color', "#7E2F8E", 'LineWidth', 1.5, 'LineStyle', '-');
    plot(y_nuc_sp, -x_nuc_sp, 'Color', "#7E2F8E", 'LineWidth', 1.5, 'LineStyle', '-');

    % Set axis properties and labels
    axis equal;
    xlim([-scalingCellLength-1.5, scalingCellLength+1.5]);
    ylim([-scalingCellWidth-1, scalingCellWidth+1]);
    title('2D Cell and Nucleoid Boundary');
    xticks([-scalingCellLength, 0, scalingCellLength]);
    xticklabels({'0','0.5','1'});
    yticks([-scalingCellWidth, 0, scalingCellWidth]);
    yticklabels({'0','0.5','1'});
    
    hold off;
end

%% -------------------------------------------------------------------------
function x_val = set_cellx(R, L, y)
% SET_CELLX Computes the x coordinate of a cell boundary for a given y.
%
%   x_val = set_cellx(R, L, y)
%
%   The cell boundary is defined piecewise with three regions:
%     - Region 1 (arc on the negative side): for y in [-(R+L), -L)
%     - Region 2 (central region): for y in [-L, L], where x is constant (R)
%     - Region 3 (arc on the positive side): for y in (L, R+L]
%
%   For y outside these regions, x_val is set to NaN.

    % Preallocate the output array
    x_val = zeros(size(y));

    % Define logical indices for each region
    region1 = (y >= -(R + L)) & (y < -L);
    region2 = (y >= -L) & (y <= L);
    region3 = (y > L) & (y <= (L + R));

    % Compute the x value for each region, ensuring non-negative arguments for sqrt.
    x_val(region1) = sqrt( max(0, R^2 - (y(region1) + L).^2) );
    x_val(region2) = R;
    x_val(region3) = sqrt( max(0, R^2 - (y(region3) - L).^2) );

    % For y outside the defined regions, assign NaN.
    x_val(~(region1 | region2 | region3)) = NaN;
end

%% -------------------------------------------------------------------------
function interpolatedVal = piecefunc(anchor, query)
% PIECEFUNC Performs piecewise linear interpolation using the provided anchor map.
%
%   interpolatedVal = piecefunc(anchor, query)
%
%   This function uses MATLAB's interp1 to compute the interpolated values

    if size(anchor, 2) < 2
        error('Anchor input must have at least two columns.');
    end

    % Perform linear interpolation (with extrapolation) using interp1.
    interpolatedVal = interp1(anchor(:,1), anchor(:,2), query, 'linear', 'extrap');
end