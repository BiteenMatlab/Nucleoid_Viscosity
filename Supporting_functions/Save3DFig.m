% Save 3D visualizations to different formats without showing figures in live script. 
% (SaveFigure crashes Matlab when saving 3D visualization to svg file)
function Save3DFig(result_folder, fileName, figHandle, fileFormats, clr_map)
    hOrig = figHandle; % Handle to the displayed figure
    hCopy = copyFigure(hOrig); % Duplicate the figure (see helper function below)
    savePlotFormats(result_folder,fileName, hCopy, ...
        fileFormats,clr_map);
end

function hCopyFig = copyFigure(hOrigFig)
% copyFigure creates a duplicate of the figure referenced by hOrigFig.
% It copies the current axes (and its children) to a new figure.
%
% The new figure is created with 'Visible' set to 'off' so that when passed to
% savePlotFormats (which will close it) the original figure remains unchanged.

% Create a new figure (hidden)
hCopyFig = figure('Visible', 'off'); 

% Get all axes from the original figure
axesHandles = findall(hOrigFig, 'Type', 'axes');


% Copy each axes to the new figure
for k = 1:length(axesHandles)
    % Copy the axes (including children) to the new figure. 
    % Note: copyobj automatically duplicates all children.
    copyobj(axesHandles(k), hCopyFig);
end

% Find the colorbar in the original figure (if any)
hColorbar = findobj(hOrigFig, 'Tag', 'Colorbar');

% If a colorbar exists, copy both the main axes and the colorbar in one call:
if ~isempty(hColorbar)
    copyobj(hColorbar, hCopyFig);
end

end

