<<<<<<< HEAD
%==========================================================================
% OrgResults.m
% -------------------------------------------------------------------------
% Consolidate workspace variables into structures containing results from
% different sections and save to Results.mat in result_folder
%==========================================================================

%% ------------------------------------------------------------------------
% 1. Localization fitting 
%% ------------------------------------------------------------------------
LocFit.PixelSize   = pixel_size;    % nm per pixel
LocFit.Conf95      = locFitCI;     % 95 % CI
LocFit.CI_UB       = locFitCI_UB;      % upper bound
clear pixel_size locFitCI locFitCI_UB

%% ------------------------------------------------------------------------
% 2. Cell Morphology & Filtering
%% ------------------------------------------------------------------------
filterValues = cell_morph(:, sel_filter);
passFilter   = filterValues > low_bd & filterValues < up_bd;
mean_ar = mean(cell_morph(passFilter, 4));    % Aspect ratio is stored in column 4.
CellMorph.AspectRatio       = mean_ar;        % aspect ratio
CellMorph.avgWidth          = cell_avg_width;
CellMorph.MorphByMovie      = morph_by_mov;
CellMorph.Filter            = sel_filter;     % morphology selected for filtering
CellMorph.filterBounds      = [low_bd, up_bd];
CellMorph.ValidIdx          = cellIndices;    % cells that pass filter
clear cell_avg_width morph_by_mov cellIndices low_bd up_bd sel_filter
clear filterValues passFilter mean_ar

%% ------------------------------------------------------------------------
% 3. Heat Maps & Spatial Displacement Maps
%% ------------------------------------------------------------------------
Maps.CellAR              = aspt_ratio;        
Maps.numYGrid            = num_pix_y;         % # pixels in y-dimension
Maps.HeatMap             = htMap;             % heat map
Maps.DispMap             = spMap;             % spatial displacement
Maps.DispCount           = spMap_cnt;         % number of steps at each grid
Maps.RawHeat             = raw_heatmap;
Maps.RawDisp             = raw_spatialDisp;
Maps.CellOutline         = cell_outline;
Maps.NucOutlineHeat      = nuc_outline;
Maps.NucOutlineDisp      = nuc_contour;
Maps.NucCntrBounds       = [cntr_low, cntr_up];
Maps.SimDispMap          = Simu_spMap;
Maps.SimSpatDisp         = Simu_spatialDisp;
clear num_pix_y htMap spMap raw_heatmap raw_spatialDisp ...
      cell_outline nuc_outline nuc_contour cntr_low cntr_up ...
      Simu_spMap Simu_spatialDisp aspt_ratio spMap_cnt

%% ------------------------------------------------------------------------
% 4. Geometry Reconstruction
%% ------------------------------------------------------------------------
% nucleoid ratio relative to cell: length, width, area
Geometry.NucToCellRatio    = Nuc_ratio_spMap;   
% cell and nucleoid parameter: length, width (normalized by cell width)
Geometry.CellParams        = cell_para;
Geometry.NucParams         = nuc_para;

Geometry.R0                = Ro;              % radius 
Geometry.L0                = Lo;              % half length of cylinder region
Geometry.ConvertFactor     = convert_factor;  % Ro/Radiu from fitting

Geometry.WidthScale        = scaling_width;   % account for cell boundary (width)
Geometry.LengthScale       = scaling_length;  % account for cell boundary (length)
Geometry.ScaleBar          = scale_bar;       % scale bar for simulation

Geometry.AnchorHeat        = anchor_htMap;    % anchor dots of nucleoid from heatmap
Geometry.NucKB_Heat        = kb_nuc_htMap;
Geometry.AnchorDisp        = anchor_spMap;    % anchor dots of nucleoid from spMap
Geometry.NucKB_Disp        = kb_nuc_spMap;

clear Nuc_ratio_spMap cell_para nuc_para Ro Lo convert_factor ...
      scaling_width scaling_length scale_bar ...
      anchor_htMap kb_nuc_htMap anchor_spMap kb_nuc_spMap

%% ------------------------------------------------------------------------
% 5. Volume & Accessibility
%% ------------------------------------------------------------------------
VolAccess.mapBounds     = [lb, rb];
VolAccess.accRange      = acc_range;
VolAccess.SimuHeatmap   = cmp_acc_htMap;
VolAccess.NucVolRatio   = nuc_volume_ratio;
VolAccess.Accessibility = acc_fit;
VolAccess.Similarity    = Metric_optz;      % Column 1: L2 norm distance (lower is better)
                                            % Column 2: Pearson correlation coefficient (higher is better)
                                            % Column 3: SSIM (higher is better)
clear lb rb acc_range cmp_acc_htMap nuc_volume_ratio acc_fit Metric_optz

%% ------------------------------------------------------------------------
% 6. Brownian-Dynamics Simulation
%% ------------------------------------------------------------------------
BDSim.NumRuns        = Num_run;
BDSim.Step           = StepSize_grid;
BDSim.Gap            = Gap_grid;
BDSim.WeightScore    = [score_w1, score_w2];
BDSim.angleAugment   = num_agl_seg;

BDSim.CoordInNuc     = coord_in_nuc;
BDSim.CoordInCyto    = coord_in_cyto;
BDSim.CoordOutCell   = coord_out_cell;

BDSim.WeightMat1     = weight_matrix1;
BDSim.WeightMat2     = weight_matrix2;

BDSim.Similarity     = Store_similarity;
BDSim.TrajFinal      = Simu_Traj_final;
BDSim.SelStepGap     = [sel_stepsize, sel_gap];

clear Num_run StepSize_grid Gap_grid score_w1 score_w2 num_agl_seg ...
      coord_in_nuc coord_in_cyto coord_out_cell ...
      weight_matrix1 weight_matrix2 Store_similarity Simu_Traj_final ...
      sel_stepsize sel_gap

%% ------------------------------------------------------------------------
% 7. Save Strucures above
%% ------------------------------------------------------------------------
save(fullfile(result_folder,'Results.mat'), ...
    'LocFit', ...
    'CellMorph', ...
    'Maps', ...
    'Geometry', ...
    'VolAccess', ...
    'BDSim', ...
    '-v7.3');             % v7.3 for large file

fprintf('\nAll structures saved to %s\n', fullfile(result_folder,'Results.mat'));

%% ------------------------------------------------------------------------
% 8. Save Raw Data (Masks, Fits, Tracks)
%% ------------------------------------------------------------------------
RawData.Mask   = Store_mask;       % binary or labeled mask(s)
RawData.Fit    = Store_fits;       % curve-fit / localization data
RawData.Track  = Store_track;      % trajectories
clear Store_mask Store_fits Store_track

save(fullfile(result_folder,'RawData.mat'), ...
    'RawData', ...
    '-v7.3');             % v7.3 for large file

fprintf('\nRaw Data saved to %s\n', fullfile(result_folder,'RawData.mat'));

%% ------------------------------------------------------------------------
% 8. Clear stored structures
%% ------------------------------------------------------------------------
=======
%==========================================================================
% OrgResults.m
% -------------------------------------------------------------------------
% Consolidate workspace variables into structures containing results from
% different sections and save to Results.mat in result_folder
%==========================================================================

%% ------------------------------------------------------------------------
% 1. Localization fitting 
%% ------------------------------------------------------------------------
LocFit.PixelSize   = pixel_size;    % nm per pixel
LocFit.Conf95      = locFitCI;     % 95 % CI
LocFit.CI_UB       = locFitCI_UB;      % upper bound
clear pixel_size locFitCI locFitCI_UB

%% ------------------------------------------------------------------------
% 2. Cell Morphology & Filtering
%% ------------------------------------------------------------------------
filterValues = cell_morph(:, sel_filter);
passFilter   = filterValues > low_bd & filterValues < up_bd;
mean_ar = mean(cell_morph(passFilter, 4));    % Aspect ratio is stored in column 4.
CellMorph.AspectRatio       = mean_ar;        % aspect ratio
CellMorph.avgWidth          = cell_avg_width;
CellMorph.MorphByMovie      = morph_by_mov;
CellMorph.Filter            = sel_filter;     % morphology selected for filtering
CellMorph.filterBounds      = [low_bd, up_bd];
CellMorph.ValidIdx          = cellIndices;    % cells that pass filter
clear cell_avg_width morph_by_mov cellIndices low_bd up_bd sel_filter
clear filterValues passFilter mean_ar

%% ------------------------------------------------------------------------
% 3. Heat Maps & Spatial Displacement Maps
%% ------------------------------------------------------------------------
Maps.CellAR              = aspt_ratio;        
Maps.numYGrid            = num_pix_y;         % # pixels in y-dimension
Maps.HeatMap             = htMap;             % heat map
Maps.DispMap             = spMap;             % spatial displacement
Maps.DispCount           = spMap_cnt;         % number of steps at each grid
Maps.RawHeat             = raw_heatmap;
Maps.RawDisp             = raw_spatialDisp;
Maps.CellOutline         = cell_outline;
Maps.NucOutlineHeat      = nuc_outline;
Maps.NucOutlineDisp      = nuc_contour;
Maps.NucCntrBounds       = [cntr_low, cntr_up];
Maps.SimDispMap          = Simu_spMap;
Maps.SimSpatDisp         = Simu_spatialDisp;
clear num_pix_y htMap spMap raw_heatmap raw_spatialDisp ...
      cell_outline nuc_outline nuc_contour cntr_low cntr_up ...
      Simu_spMap Simu_spatialDisp aspt_ratio spMap_cnt

%% ------------------------------------------------------------------------
% 4. Geometry Reconstruction
%% ------------------------------------------------------------------------
% nucleoid ratio relative to cell: length, width, area
Geometry.NucToCellRatio    = Nuc_ratio_spMap;   
% cell and nucleoid parameter: length, width (normalized by cell width)
Geometry.CellParams        = cell_para;
Geometry.NucParams         = nuc_para;

Geometry.R0                = Ro;              % radius 
Geometry.L0                = Lo;              % half length of cylinder region
Geometry.ConvertFactor     = convert_factor;  % Ro/Radiu from fitting

Geometry.WidthScale        = scaling_width;   % account for cell boundary (width)
Geometry.LengthScale       = scaling_length;  % account for cell boundary (length)
Geometry.ScaleBar          = scale_bar;       % scale bar for simulation

Geometry.AnchorHeat        = anchor_htMap;    % anchor dots of nucleoid from heatmap
Geometry.NucKB_Heat        = kb_nuc_htMap;
Geometry.AnchorDisp        = anchor_spMap;    % anchor dots of nucleoid from spMap
Geometry.NucKB_Disp        = kb_nuc_spMap;

clear Nuc_ratio_spMap cell_para nuc_para Ro Lo convert_factor ...
      scaling_width scaling_length scale_bar ...
      anchor_htMap kb_nuc_htMap anchor_spMap kb_nuc_spMap

%% ------------------------------------------------------------------------
% 5. Volume & Accessibility
%% ------------------------------------------------------------------------
VolAccess.mapBounds     = [lb, rb];
VolAccess.accRange      = acc_range;
VolAccess.SimuHeatmap   = cmp_acc_htMap;
VolAccess.NucVolRatio   = nuc_volume_ratio;
VolAccess.Accessibility = acc_fit;
VolAccess.Similarity    = Metric_optz;      % Column 1: L2 norm distance (lower is better)
                                            % Column 2: Pearson correlation coefficient (higher is better)
                                            % Column 3: SSIM (higher is better)
clear lb rb acc_range cmp_acc_htMap nuc_volume_ratio acc_fit Metric_optz

%% ------------------------------------------------------------------------
% 6. Brownian-Dynamics Simulation
%% ------------------------------------------------------------------------
BDSim.NumRuns        = Num_run;
BDSim.Step           = StepSize_grid;
BDSim.Gap            = Gap_grid;
BDSim.WeightScore    = [score_w1, score_w2];
BDSim.angleAugment   = num_agl_seg;

BDSim.CoordInNuc     = coord_in_nuc;
BDSim.CoordInCyto    = coord_in_cyto;
BDSim.CoordOutCell   = coord_out_cell;

BDSim.WeightMat1     = weight_matrix1;
BDSim.WeightMat2     = weight_matrix2;

BDSim.Similarity     = Store_similarity;
BDSim.TrajFinal      = Simu_Traj_final;
BDSim.SelStepGap     = [sel_stepsize, sel_gap];

clear Num_run StepSize_grid Gap_grid score_w1 score_w2 num_agl_seg ...
      coord_in_nuc coord_in_cyto coord_out_cell ...
      weight_matrix1 weight_matrix2 Store_similarity Simu_Traj_final ...
      sel_stepsize sel_gap

%% ------------------------------------------------------------------------
% 7. Save Strucures above
%% ------------------------------------------------------------------------
save(fullfile(result_folder,'Results.mat'), ...
    'LocFit', ...
    'CellMorph', ...
    'Maps', ...
    'Geometry', ...
    'VolAccess', ...
    'BDSim', ...
    '-v7.3');             % v7.3 for large file

fprintf('\nAll structures saved to %s\n', fullfile(result_folder,'Results.mat'));

%% ------------------------------------------------------------------------
% 8. Save Raw Data (Masks, Fits, Tracks)
%% ------------------------------------------------------------------------
RawData.Mask   = Store_mask;       % binary or labeled mask(s)
RawData.Fit    = Store_fits;       % curve-fit / localization data
RawData.Track  = Store_track;      % trajectories
clear Store_mask Store_fits Store_track

save(fullfile(result_folder,'RawData.mat'), ...
    'RawData', ...
    '-v7.3');             % v7.3 for large file

fprintf('\nRaw Data saved to %s\n', fullfile(result_folder,'RawData.mat'));

%% ------------------------------------------------------------------------
% 8. Clear stored structures
%% ------------------------------------------------------------------------
>>>>>>> 5e4cf7b (Initial commit)
clear LocFit CellMorph Maps Geometry VolAccess BDSim RawData