% % avoid showing figure in live script after saving figures
function SaveFigure(result_folder, fileName, figHandle, fileFormats, clr_map)
% SaveFigure  Supports PNG, JPG, TIFF, PDF, EPS, SVG, EMF, FIG.

% ------------------------------------------------------------------
% 1. Clone the entire figure (keeps colour-bars, legends, layout)
% ------------------------------------------------------------------
hCopy = copyobj(figHandle, 0);          % duplicate under the root
set(hCopy,'Visible','off');

% Keep same on-screen pixel size
origUnits = get(figHandle,'Units');      set(figHandle,'Units','pixels');
pixPos    = get(figHandle,'Position');   set(hCopy,'Units','pixels','Position',pixPos);
set(figHandle,'Units',origUnits);        set(hCopy,'PaperPositionMode','auto');

% Optional colormap override
if nargin>=5 && ~isempty(clr_map)
    colormap(hCopy, clr_map);
end

% ------------ critical layout pass (fixes overlaps) ----------------
set(0,'CurrentFigure',hCopy);   % make duplicate current
drawnow update;                 % finish tiledlayout/sgtitle adjustments

% ------------------------------------------------------------------
% 2. Export
% ------------------------------------------------------------------
for k = 1:numel(fileFormats)
    fmt   = lower(fileFormats{k});
    fname = fullfile(result_folder, [fileName '.' fmt]);

    switch fmt
        case {'png','jpg','jpeg','tif','tiff','bmp'}
            exportgraphics(hCopy, fname);

        case 'pdf'
            exportgraphics(hCopy, fname, 'ContentType','vector');

        case 'eps'
            print(hCopy, fname, '-depsc', '-vector');

        case 'svg'
            try
                exportgraphics(hCopy, fname, 'ContentType','vector');
            catch
                print(hCopy, fname, '-dsvg', '-vector');
            end

        case 'emf'                      % Windows only
            print(hCopy, fname, '-dmeta', '-vector');

        case 'fig'
            savefig(hCopy, fname);

        otherwise
            warning('SaveFigure:UnknownFormat', ...
                    'Unknown format "%s" — skipping.', fmt);
    end
end

close(hCopy);
end

