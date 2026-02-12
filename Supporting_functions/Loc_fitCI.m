function [locPrecision, locPrecisionUB] = Loc_fitCI(Store_track, Store_fits, pixel_size, outlier_method)
% Loc_fitCI Calculates the localization fitting confidence interval of single-particle tracking data.
%
% [locPrecision, locPrecisionUB] = Loc_fitCI(Store_track, Store_fits, pixel_size, outlier_method)
%
% This function extracts confidence interval (CI) values from the fit data (Store_fits)
% corresponding to the molecule indices in the tracking data (Store_track). Following
% concatenation of all CI values (in pixels), outliers are removed using the specified
% outlier detection method (default is 'mean'). A diagnostic histogram is plotted.
%
% The mean CI is used as the localization precision and is scaled by pixel_size to convert
% to physical units. Similarly, the maximum CI is computed as an upper bound.
%
% Inputs:
% Store_track - Cell array where each element is an [N x M] numeric matrix containing
% tracking data. The 6th column should contain molecule indices (positive integers).
% Store_fits - Cell array (same length as Store_track) where each element is a structure
% with a numeric field called 'widthcCI' containing confidence interval values.
% pixel_size - Positive scalar used to convert CI values from pixels to physical units.
% outlier_method- (Optional) String or method accepted by isoutlier (e.g., 'mean', 'median', etc.).
% Default is 'mean'.
%
% Outputs:
% locPrecision - Mean localization precision in physical units.
% locPrecisionUB - Maximum localization precision (upper bound) in physical units.
%
% Example:
% [lp, lpUB] = Loc_fitCI(Store_track, Store_fits, 49, 'median');
%
% Author: Xiaofeng Dai
% Date: 05/05/2025

%% 1. Input Validation and Default Setup
if nargin < 4 || isempty(outlier_method)
    outlier_method = 'mean';
end
if ~iscell(Store_track) || ~iscell(Store_fits)
    error('Store_track and Store_fits must be cell arrays.');
end
if numel(Store_track) ~= numel(Store_fits)
    error('Store_track and Store_fits must have the same number of elements.');
end
if ~isscalar(pixel_size) || pixel_size <= 0
    error('pixel_size must be a positive scalar.');
end

numFiles = numel(Store_track);
ciValues = cell(numFiles, 1);

%% 2. Parallel Extraction of Confidence Interval Values
parfor i = 1:numFiles
    currTrack = Store_track{i};
    currFit = Store_fits{i};

    % Validate the tracking data structure; must have at least 6 columns.
    if size(currTrack, 2) < 6
        error('Element %d of Store_track must have at least 6 columns.', i);
    end

    % Ensure the fit structure contains the field "widthcCI".
    if ~isfield(currFit, 'widthcCI')
        error('Store_fits{%d} does not contain the field "widthcCI".', i);
    end

    % Extract the molecule indices (col 6 from SMALLLABS analysis output).
    moleculeIndices = currTrack(:,6);
    widthCI = currFit.widthcCI;

    % Verify that the molecule indices are valid
    if any(moleculeIndices < 1) || any(moleculeIndices > numel(widthCI))
        error('Invalid molecule index in Store_track{%d}.', i);
    end

    % Extract the CI values using molecule indices.
    ciValues{i} = widthCI(moleculeIndices);
end

%% 3. Concatenate and Clean CI Values
allCI = cat(1, ciValues{:});
% Remove outliers according to the specified method.
allCI(isoutlier(allCI, outlier_method)) = [];
% Remove any NaN values.
allCI(isnan(allCI)) = [];

%% 4. Plot Diagnostic Histogram
figure;
histogram(allCI);
title('Distribution of Localization Confidence Intervals');
xlabel('Width Confidence Interval (pixels)');
ylabel('Probability');

%% 5. Compute Localization Precision (in physical units)
locPrecision   = mean(allCI) * pixel_size;
locPrecisionUB = max(allCI) * pixel_size;

%% 6. Display Diagnostic Information
fprintf('Mean localization Confidence Interval = %.3f nm \n', locPrecision);