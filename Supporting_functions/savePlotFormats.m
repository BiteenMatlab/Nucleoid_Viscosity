function savePlotFormats(resultFolder, baseFileName, figHandle, fileFormats, clr_map)
% savePlotFormats saves a figure into multiple file formats.
%
% savePlotFormats(resultFolder, baseFileName) saves the current figure (gcf)
% into the folder specified by resultFolder with a file name starting with
% baseFileName in the default formats: .fig, .png, and .svg. Using 'parula'
% as default color map
%
% savePlotFormats(resultFolder, baseFileName, figHandle) saves the provided
% figure handle.
%
% savePlotFormats(resultFolder, baseFileName, figHandle, fileFormats) saves the
% figure in the file formats (extensions) given by the cell array fileFormats.
% For example, to save as .fig, .png, and .svg set fileFormats = {'fig','png','svg'}.
%
% INPUTS:
% resultFolder - (String) Folder in which to save the files. It must exist.
% baseFileName - (String) Base name for the saved figure files (without extension).
% figHandle - (Optional, Figure handle) The figure to save.
% Default: gcf.
% fileFormats - (Optional, Cell array of strings) File formats to save.
% Default: {'fig','png','svg'}.
% clr_map - (Optional) choose color map.
% Default: 'parula'.
%
% EXAMPLES:
% % Save the current figure into default formats:
% savePlotFormats('C:\Results','Cell Morohlogy Distribution');
%
% % Save a given figure into .png and .pdf, with 'hot' as color map:
% savePlotFormats('C:\Results','Cell Morohlogy Distribution', myFigure, {'png','pdf'},'hot');
%
% Author: Xiaofeng Dai
% Date: 04/08/2025

%% 1. Input Validation and Defaults
narginchk(2, 5); % At least 2 inputs required

% If figHandle is not provided or empty, use the current figure
if nargin < 3 || isempty(figHandle)
    figHandle = gcf;
end

% Set default fileFormats if not provided
if nargin < 4 || isempty(fileFormats)
    fileFormats = {'fig', 'png', 'svg'};
end

% Set default color map if not provided
if nargin < 5 || isempty(clr_map)
    clr_map = 'parula';
end

% Validate that resultFolder exists
if ~isfolder(resultFolder)
    error('The result folder "%s" does not exist.', resultFolder);
end

% Validate that baseFileName is a string/char
if ~(ischar(baseFileName) || isstring(baseFileName))
    error('baseFileName must be a character vector or a string.');
end

% Optionally adjust figure properties (set position/size)
figHandle.Position = [100, 100, 1000, 800];
colormap(figHandle,clr_map)

%% 2. Save Figure in Each Assigned Format
% Loop through the list of file formats and save accordingly.
for k = 1:numel(fileFormats)
    % Remove any leading '.' if present and force lower case.
    fmt = lower(regexprep(fileFormats{k}, '^\.', ''));
    % Build the complete file name with extension.
    fullFileName = fullfile(resultFolder, [char(baseFileName), '.', fmt]);

    switch fmt
        case 'fig'
            % Use savefig for MATLAB figure files.
            savefig(figHandle, fullFileName);
        otherwise
            % For formats such as png, svg, pdf, etc., use saveas.
            saveas(figHandle, fullFileName);
    end
end
close(figHandle);
end