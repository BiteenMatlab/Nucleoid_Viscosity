function [R, L, residual] = Rod_Morph_fit(outline, num_pix_y, aspt_ratio, figshow)
% Rod_Morph_fit Fits a rod‐shaped model (a central cylinder with two semicircular ends)
% to the upper half of a cell boundary outline.
%
% [R, L, residual] = Rod_Morph_fit(outline, num_pix_y, aspt_ratio, figshow)
%
% Inputs:
% outline - [N x 2] numeric array of (x,y) coordinates representing the cell
% boundary. The model will be fitted to the upper half only.
% num_pix_y - Scalar indicating the number of grid pixels (y–dimension) used
% in the normalization (e.g., 20).
% aspt_ratio - Scalar specifying the aspect ratio (ratio of x–dimension to y–dimension),
% e.g., 2.
% figshow - (Optional) 'on' to display a plot of the fit; 'off' (default) otherwise.
%
% Outputs:
% R - The fitted cell half-width (rod radius).
% L - The half-length of the central cylindrical part.
% residual - The mean squared residual of the fit (i.e., LSQ residual normalized by the
% number of data points).
%
% The rod model is defined piecewise for the top half:
% For x in [-(R+L), -L]: y = sqrt(R^2 - (x+L)^2)
% For x in [-L, L]: y = R
% For x in [L, R+L]: y = sqrt(R^2 - (x-L)^2)
%
% The function first shifts the outline so that the model’s center is at (0,0)
% (using provided normalization parameters) and then fits the positive (upper)
% portion only.
%
% Example:
% [R, L, err] = Rod_Morph_fit(cellOutline, 20, 2, 'on');
%
% Author: Xiaofeng Dai
% Date: 04/14/2025

%% Input Handling and Defaults
narginchk(3,4);
if nargin < 4 || isempty(figshow)
    figshow = 'off';
end
% Validate outline dimensions.
if ~ismatrix(outline) || size(outline,2) ~= 2
    error('Input "outline" must be an [N x 2] numeric array.');
end

%% 1. Shift the Outline to Center the Rod Model
% Compute shift values: the expected "center" is at (num_pix_y*aspt_ratio/2 + 0.5, num_pix_y/2+0.5).
shiftx = num_pix_y * aspt_ratio / 2 + 0.5;
shifty = num_pix_y / 2 + 0.5;
% Shift coordinates so that the center becomes [0,0].
outline(:,1) = outline(:,1) - shiftx;
outline(:,2) = outline(:,2) - shifty;
% Make a copy (for later plotting) of the shifted outline.
buff = outline;

%% 2. Use Only the Upper Half of the Outline
% (Assuming symmetry: fit the rod to points with y >= 0.)
outline = outline(outline(:,2) >= 0, :);
% If no points remain, error out.
if isempty(outline)
    error('No valid outline points with nonnegative y were found after shifting.');
end

%% 3. Define the Rod Model Function (Upper Half)
% F(Morph,x) returns the model y for each x.
% Here, Morph(1) = R (radius) and Morph(2) = L (half-length of central cylinder).
F = @(Morph, x) ((x >= -(Morph(1)+Morph(2)) & x < -Morph(2)) .* sqrt( max(0, Morph(1)^2 - (x + Morph(2)).^2) ) + ...
                 (x >= -Morph(2) & x <= Morph(2)) .* Morph(1) + ...
                 (x >  Morph(2) & x <= (Morph(1)+Morph(2))) .* sqrt( max(0, Morph(1)^2 - (x - Morph(2)).^2) ));

%% 4. Choose Initial Guess and Bounds for Fitting
% Initial guess: R from maximum y value; L roughly as (max_x - R) [assumes positive x values].
init_R = max(outline(:,2));
init_L = max(outline(:,1)) - init_R;
morph0 = [init_R, init_L];
lb = [0, 0];
ub = [num_pix_y/2, num_pix_y * aspt_ratio / 2];

%% 5. Fit the Model using lsqcurvefit
opt = optimset('Display', 'off');
[RL, resnorm] = lsqcurvefit(F, morph0, outline(:,1), outline(:,2), lb, ub, opt);
R = RL(1);
L = RL(2);
% Compute mean squared residual (number of points used).
residual = resnorm / size(outline,1);

%% 6. (Optional) Plot the Fitted Rod Model Versus Data
if strcmpi(figshow, 'on')
    % Define an inline version of the rod model for plotting.
    Outfit = @(y) ((y >= -(R+L) & y < -L) .* sqrt( max(0, R^2 - (y + L).^2) ) + ...
                   (y >= -L & y <= L) .* R + ...
                   (y > L & y <= (R+L)) .* sqrt( max(0, R^2 - (y - L).^2) ));
    % Generate test x–values spanning the rod width.
    testy = -R - L : 0.1 : R + L;
    % Suppress warnings about small negative values under the square-root.
    warning('off', 'all');
    testx = Outfit(testy);
    figure;
    % Plot upper half (blue line) and mirror for the lower half (blue dashed line).
    plot(testy, testx, 'b', 'LineWidth', 1.5);
    hold on;
    plot(testy, -testx, 'b', 'LineWidth', 1.5);
    hold on;
    % Plot the data points (shifted outline; use buff from before filtering).
    scatter(buff(:,1), buff(:,2), 10, 'red', 'filled');
    axis equal;
    xlim([-(num_pix_y*aspt_ratio/2),num_pix_y*aspt_ratio/2])
    ylim([-num_pix_y/2,num_pix_y/2])
    title('Rod Fitting');
    warning('on', 'all');

    fprintf('Rfit = %.2f;  Lfit = %.2f.\nThe mean squared error is %.5f\n', R, L, residual);
end
end