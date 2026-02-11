function vars = Results4Btstrp(result_folder)

% -------------------------------------------------------------------------
% 0. Locate folder
% -------------------------------------------------------------------------
if nargin < 1 || isempty(result_folder)
    result_folder = uigetdir(pwd,'Select folder');
    if result_folder == 0
        disp('No folder was selected.');  return
    end
end

% -------------------------------------------------------------------------
% 1. Read the Results.mat 
% -------------------------------------------------------------------------
R = load(fullfile(result_folder,'Results.mat'), ...
          'CellMorph','LocFit','Maps');

% Pull the requested fields
vars.cellIndices     = R.CellMorph.ValidIdx;
vars.morph_by_mov    = R.CellMorph.MorphByMovie;
vars.num_pix_y       = R.Maps.numYGrid;
vars.locFitCI_UB     = R.LocFit.CI_UB;
vars.pixel_size      = R.LocFit.PixelSize;
vars.cell_avg_width  = R.CellMorph.avgWidth;

% -------------------------------------------------------------------------
% 2. Read RawData.mat
% -------------------------------------------------------------------------
RD = load(fullfile(result_folder,'RawData.mat'),'RawData');
vars.Store_track = RD.RawData.Track;
vars.Store_mask  = RD.RawData.Mask;
vars.Store_fits  = RD.RawData.Fit;

% -------------------------------------------------------------------------
% 3. Either return the struct or splash into workspace
% -------------------------------------------------------------------------
if nargout == 0   
    wsc = 'caller';
    names = fieldnames(vars);
    for k = 1:numel(names)
        assignin(wsc, names{k}, vars.(names{k}));
    end
    clearvars vars
else
    % just return the struct
end
end
