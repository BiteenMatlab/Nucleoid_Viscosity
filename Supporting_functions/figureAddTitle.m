<<<<<<< HEAD
function figureAddTitle(figHandle, titlename, low_bd, up_bd, aspt_ratio)
% figureAddTitle Adds a composite title to a figure.
%
% figureAddTitle(figHandle, low_bd, up_bd, aspt_ratio) adds a main title and
% subtitle to the current axes of the specified figure. The subtitles include
% the lower bound, upper bound, and aspect ratio information.
%
% Inputs:
% figHandle - Handle to a figure in which the title will be added.
% If empty or not provided, gcf is used.
% titlename - title of the plot.
% low_bd - Numeric lower bound (e.g., for filtering), which will be shown
% in the subtitle.
% up_bd - Numeric upper bound.
% aspt_ratio- Numeric aspect ratio.
%
% The function creates a title with:
% Main Title: 'Localization Heat Map after filtering'
% Subtitle: ' Lower bound = µm, Upper bound = µm. Aspect ratio = '
%
% It then adjusts the vertical positions of the title and subtitle.
%
% Example:
% figureAddTitle(gcf, 1.2, 2.2, 2);
%
% Author: Xiaofeng Dai
% Date: 04/12/2025

%% 1. Validate Inputs and Get Current Axes
if nargin < 1 || isempty(figHandle)
    figHandle = gcf;
end
% Get the current axes; create one if none exists.
ax = get(figHandle, 'CurrentAxes');
if isempty(ax)
    ax = axes(figHandle);
end

%% 2. Format the Subtitle String
% Create individual strings for lower bound, upper bound, and aspect ratio.
% (A leading space is added for cosmetic reasons.)
str1 = sprintf(' Lower bound = %.1f', low_bd);
str2 = sprintf('Upper bound = %.1f', up_bd);
str3 = sprintf('Aspect ratio = %.2f', aspt_ratio);
% Combine into one subtitle string.
annotationStr = [str1 ' µm, ' str2 ' µm. ' str3];

%% 3. Set the Title with Main Title and Subtitle
[t, s] = title(ax, titlename, annotationStr);

%% 4. Adjust the Title Positions
% Retrieve the positions for the title objects.
posMain = get(t, 'position');  % [x, y, z] of main title
posSub  = get(s, 'position');   % [x, y, z] of subtitle
% Adjust the vertical position (y component) to fine-tune the placements.
posMain = posMain + [0, 0.27, 0];
posSub  = posSub  + [0, 0.145, 0];
set(t, 'position', posMain);
set(s, 'position', posSub);
=======
function figureAddTitle(figHandle, titlename, low_bd, up_bd, aspt_ratio)
% figureAddTitle Adds a composite title to a figure.
%
% figureAddTitle(figHandle, low_bd, up_bd, aspt_ratio) adds a main title and
% subtitle to the current axes of the specified figure. The subtitles include
% the lower bound, upper bound, and aspect ratio information.
%
% Inputs:
% figHandle - Handle to a figure in which the title will be added.
% If empty or not provided, gcf is used.
% titlename - title of the plot.
% low_bd - Numeric lower bound (e.g., for filtering), which will be shown
% in the subtitle.
% up_bd - Numeric upper bound.
% aspt_ratio- Numeric aspect ratio.
%
% The function creates a title with:
% Main Title: 'Localization Heat Map after filtering'
% Subtitle: ' Lower bound = µm, Upper bound = µm. Aspect ratio = '
%
% It then adjusts the vertical positions of the title and subtitle.
%
% Example:
% figureAddTitle(gcf, 1.2, 2.2, 2);
%
% Author: Xiaofeng Dai
% Date: 04/12/2025

%% 1. Validate Inputs and Get Current Axes
if nargin < 1 || isempty(figHandle)
    figHandle = gcf;
end
% Get the current axes; create one if none exists.
ax = get(figHandle, 'CurrentAxes');
if isempty(ax)
    ax = axes(figHandle);
end

%% 2. Format the Subtitle String
% Create individual strings for lower bound, upper bound, and aspect ratio.
% (A leading space is added for cosmetic reasons.)
str1 = sprintf(' Lower bound = %.1f', low_bd);
str2 = sprintf('Upper bound = %.1f', up_bd);
str3 = sprintf('Aspect ratio = %.2f', aspt_ratio);
% Combine into one subtitle string.
annotationStr = [str1 ' µm, ' str2 ' µm. ' str3];

%% 3. Set the Title with Main Title and Subtitle
[t, s] = title(ax, titlename, annotationStr);

%% 4. Adjust the Title Positions
% Retrieve the positions for the title objects.
posMain = get(t, 'position');  % [x, y, z] of main title
posSub  = get(s, 'position');   % [x, y, z] of subtitle
% Adjust the vertical position (y component) to fine-tune the placements.
posMain = posMain + [0, 0.27, 0];
posSub  = posSub  + [0, 0.145, 0];
set(t, 'position', posMain);
set(s, 'position', posSub);
>>>>>>> 5e4cf7b (Initial commit)
end