function Geometry_Reconstruct_3D(Radius, Length, scalingCellWidth, scalingCellLength, anchor_nuc)
% GEOMETRY_RECONSTRUCT_3D Visualizes the 3D geometry of reconstructed cell and nucleoid.
%
%   Geometry_Reconstruct_3D(Radius, Length, scalingCellWidth, scalingCellLength, anchor_nuc)
%
%   Inputs:
%     Radius            - Effective radius for the spherical caps and cylinder.
%     Length            - Length offset for the rod-shaped cell (distance from midcell to cap centers).
%     scalingCellWidth  - Scale for the cell’s width (used to set x and z axis limits).
%     scalingCellLength - Scale for the cell’s length (used to set y axis limits).
%     anchor_nuc        - An [n x 2] anchor map for the nucleoid boundary.
%
%   This function computes the 3D cell shape using a rod model (spherical end-caps plus
%   a cylindrical middle section) and the nucleoid shape using a rotated profile. 
%   The results are displayed in three subplots.
%
% Author: Xiaofeng Dai
% Date: 05/11/2025

    %% ------- Input Validation -------
    narginchk(5,5);
    
    if size(anchor_nuc,2) < 2
        error('anchor_nuc must have at least two columns for coordinates.');
    end
    if size(anchor_nuc,1) < 3
        error('anchor_nuc must have at least three anchor points.');
    end
    if ~issorted(anchor_nuc(:,1))
        warning('anchor_nuc first column is not sorted. Sorting anchor_nuc now.');
        anchor_nuc = sortrows(anchor_nuc,1);
    end

    %% ------- Define Nucleoid Boundary Function -------
    % Function handle that returns nucleoid boundary based on anchor_nuc.
    nucBoundaryFun = @(y) piecefunc(anchor_nuc, y);

    %% ------- Compute 3D Cell Shape -------
    % Get the 3D coordinates for the rod-shaped cell.
    [X_cell, Y_cell, Z_cell] = Rod3D(Radius, Length);
    % Combine the three components into a single 3D array.
    surf_cell = cat(3, X_cell, Y_cell, Z_cell);
    
    %% ------- Visualization -------
    figure;
    figure('Units','pixels','Position',[100 100 900 2200]);
    % --- Subplot 1: Cell Shape ---
    subplot(3,1,1);
    mesh(surf_cell(:,:,1), surf_cell(:,:,2), surf_cell(:,:,3));
    axis equal;
    title('Cell Shape');
    xlim([-scalingCellWidth-1, scalingCellWidth+1]);
    ylim([-scalingCellLength-1.5, scalingCellLength+1.5]);
    zlim([-scalingCellWidth-1, scalingCellWidth+1]);
    xticks([-scalingCellWidth, 0, scalingCellWidth]);
    xticklabels({'-0.5','0','0.5'});
    % yticks([-scalingCellLength, 0, scalingCellLength]);
    % yticklabels({'-0.5','0','0.5'});
    yticks([-scalingCellLength,-0.5*scalingCellLength, 0,0.5*scalingCellLength, scalingCellLength]);
    yticklabels({'-0.5','-0.25','0','0.25','0.5'});
    zticks([-scalingCellWidth, 0, scalingCellWidth]);
    zticklabels({'-0.5','0','0.5'});
    
    % --- Subplot 2: Nucleoid Shape ---
    % Generate the nucleoid shell.
    surf_nuc = Shell3D(nucBoundaryFun, anchor_nuc);
    subplot(3,1,2);
    mesh(surf_nuc(:,:,1), surf_nuc(:,:,2), surf_nuc(:,:,3), 'EdgeColor', 'interp');
    axis equal;
    title('Nucleoid Shape');
    
    % Set axis limits based on anchor_nuc extents.
    l_end = max(anchor_nuc(:,1));
    w_end = max(anchor_nuc(:,2));
    if l_end - floor(l_end) > 0.75
        y_tick = ceil(l_end);
    else
        y_tick = floor(l_end);
    end
    if w_end - floor(w_end) > 0.75
        x_tick = ceil(w_end);
    else
        x_tick = floor(w_end);
    end
    ylim([-l_end-0.5, l_end+0.5]);
    xlim([-w_end-0.25, w_end+0.25]);
    zlim([-w_end-0.25, w_end+0.25]);
    yticks([-y_tick, 0, y_tick]);
    y_label = round(y_tick / scalingCellLength / 2, 2);
    yticklabels({string(-y_label), '0', string(y_label)});
    xticks([-x_tick, 0, x_tick]);
    x_label = round(x_tick / scalingCellWidth / 2, 2);
    xticklabels({string(-x_label), '0', string(x_label)});
    zticks([-x_tick, 0, x_tick]);
    zticklabels({string(-x_label), '0', string(x_label)});
    
    % --- Subplot 3: Combined Cell and Nucleoid ---
    % For better visualization, mask the upper half of the cell surface (where Z > 0).
    idx = find(Z_cell > 0);
    X_cell(idx) = NaN;
    Y_cell(idx) = NaN;
    Z_cell(idx) = NaN;
    
    subplot(3,1,3);
    mesh(X_cell, Y_cell, Z_cell);
    hold on;
    mesh(surf_nuc(:,:,1), surf_nuc(:,:,2), surf_nuc(:,:,3), 'EdgeColor', 'interp');
    axis equal;
    title('Cell and Nucleoid Visualization');
    xlim([-scalingCellWidth-1, scalingCellWidth+1]);
    ylim([-scalingCellLength-1.5, scalingCellLength+1.5]);
    zlim([-scalingCellWidth-1, scalingCellWidth+1]);
    xticks([-scalingCellWidth, 0, scalingCellWidth]);
    xticklabels({'-0.5','0','0.5'});
    % yticks([-scalingCellLength, 0, scalingCellLength]);
    % yticklabels({'-0.5','0','0.5'});
    yticks([-scalingCellLength,-0.5*scalingCellLength, 0,0.5*scalingCellLength, scalingCellLength]);
    yticklabels({'-0.5','-0.25','0','0.25','0.5'});
    zticks([-scalingCellWidth, 0, scalingCellWidth]);
    zticklabels({'-0.5','0','0.5'});
    hold off;
    
    % Clear temporary variables.
    clear X_cell Y_cell Z_cell idx surf_nuc surf_cell;
end

%% -------------------------------------------------------------------------
function Shell_coord = Shell3D(fun, anchor)
% SHELL3D Generates a 3D shell surface.
%
%   Shell_coord = Shell3D(fun, anchor)
%
%   Inputs:
%     fun   - Function handle that computes the radial boundary given y values.
%     anchor- An [n x 2] anchor map; the first column contains y values.
%
%   Output:
%     Shell_coord - A 3D array of size [numTheta x numY x 3] with the X, Y, Z
%                   coordinates of the nucleoid shell surface.
%
%   The radial profile (computed via fun) is rotated around the y-axis to produce a shell.
%

    % Define y values based on the anchor map.
    y_range = linspace(min(anchor(:,1)), max(anchor(:,1)), 500);
    
    % Compute the radial distances for these y values.
    radial_positive = double(fun(y_range));
    radial_negative = -radial_positive;
    
    % Combine the positive and negative radial profiles.
    radial_combined = [radial_positive; radial_negative];  % [2 x length(y_range)]
    
    % Define the angular coordinate (theta) for revolution.
    theta_vals = linspace(0, 2*pi, 200);

    Y_mesh = cat(1,y_range,y_range);
    Y_mesh = meshgrid(Y_mesh,theta_vals);

    [X_rot,theta] = meshgrid(radial_combined,theta_vals);
    
    % Compute the 3D coordinates via rotation about the y-axis.
    X = X_rot.*cos(theta);
    Z = X_rot.*sin(theta);
    
    
    % Assemble the coordinates into a 3D array.
    Shell_coord = cat(3, X, Y_mesh, Z);
    Shell_coord = real(Shell_coord);
end

%% -------------------------------------------------------------------------
function [X_all, Y_all, Z_all] = Rod3D(R, L)
% ROD3D Constructs a 3D rod-shaped geometry representing a cell.
%
%   [X_all, Y_all, Z_all] = Rod3D(R, L)
%
%   Inputs:
%     R - Radius of the spherical caps and cylinder.
%     L - Length offset for the rod (distance between the midplane and cap centers).
%
%   Outputs:
%     X_all, Y_all, Z_all - Matrices representing the 3D coordinates of the rod surface.
%
%   The rod is constructed by combining:
%     1. The left spherical cap,
%     2. A middle cylindrical section,
%     3. The right spherical cap.
%

    % Create a grid for generating spherical caps.
    gridPts = 501;
    x_vals = linspace(-R, R, gridPts);
    z_vals = linspace(-R, R, gridPts);
    [X, Z] = meshgrid(x_vals, z_vals);
    
    % Compute the upper spherical cap (shifted upward by L).
    % Use max(0,·) to ensure a non-negative argument inside sqrt.
    Y_left = sqrt(max(0, R^2 - X.^2 - Z.^2)) + L;
    % Mask points where the argument is negative.
    Y_left(R^2 - X.^2 - Z.^2 < 0) = NaN;
    
    % Compute the lower spherical cap (shifted downward by L).
    Y_right = -sqrt(max(0, R^2 - X.^2 - Z.^2)) - L;
    Y_right(R^2 - X.^2 - Z.^2 < 0) = NaN;
    
    % Create the cylindrical section.
    [X_cyl, Z_cyl, Y_cyl] = cylinder(R, length(X)-1);
    % Adjust cylinder height to span the gap between the two caps.
    Y_cyl = (Y_cyl - 0.5) * (2 * L+1);
    
    % Concatenate the three parts.
    % It is assumed that the X and Z grids for the spherical caps are compatible.
    X_all = [X; X_cyl; X];
    Y_all = [Y_left; Y_cyl; Y_right];
    Z_all = [Z; Z_cyl; Z];
end

%% -------------------------------------------------------------------------
function interpolatedVal = piecefunc(anchor, query)
% PIECEFUNC Performs piecewise linear interpolation using the provided anchor map.
%
%   interpolatedVal = piecefunc(anchor, query)
%
%   The anchor input must be a two-column matrix [independent, dependent]. This function uses
%   MATLAB's interp1 to compute the interpolated (or extrapolated) values at the query points.
%

    if size(anchor, 2) < 2
        error('Anchor input must have at least two columns [independent, dependent].');
    end
    
    % Compute the interpolated value with linear interpolation and extrapolation.
    interpolatedVal = interp1(anchor(:,1), anchor(:,2), query, 'linear', 'extrap');
end