<<<<<<< HEAD
function UnpackResults(result_folder, rawData)

% -------------------------------------------------------------------------
% 0. Locate results --------------------------------------------------------
% -------------------------------------------------------------------------
if nargin < 1 || isempty(result_folder)
    result_folder = uigetdir(pwd, 'Select result folder');
    if result_folder == 0
        disp('No folder was selected');
        return
    end
end

if nargin < 2 || isempty(rawData)
    rawData = 'off';
end


% -------------------------------------------------------------------------
% 1. Load structures -------------------------------------------------------
% -------------------------------------------------------------------------
load(fullfile(result_folder, 'Results.mat'), ...
     'LocFit','CellMorph','Maps','Geometry','VolAccess','BDSim');

% Put everything into the caller workspace so the variables behave as if
% they had just been computed there.
wsc = 'caller';

% -------------------------------------------------------------------------
% 2. LocFit ---------------------------------------------------------------
assignin(wsc, 'pixel_size',  LocFit.PixelSize);
assignin(wsc, 'locFitCI',    LocFit.Conf95);
assignin(wsc, 'locFitCI_UB', LocFit.CI_UB);
clear LocFit

% -------------------------------------------------------------------------
% 3. CellMorph & filtering info -------------------------------------------
assignin(wsc, 'mean_ar',         CellMorph.AspectRatio);
assignin(wsc, 'cell_avg_width',  CellMorph.avgWidth);
assignin(wsc, 'morph_by_mov',    CellMorph.MorphByMovie);
assignin(wsc, 'sel_filter',      CellMorph.Filter);
assignin(wsc, 'low_bd',          CellMorph.filterBounds(1));
assignin(wsc, 'up_bd',           CellMorph.filterBounds(2));
assignin(wsc, 'cellIndices',     CellMorph.ValidIdx);
clear CellMorph

% -------------------------------------------------------------------------
% 4. Heat‑ & displacement maps --------------------------------------------
assignin(wsc, 'aspt_ratio',        Maps.CellAR);
assignin(wsc, 'num_pix_y',         Maps.numYGrid);
assignin(wsc, 'htMap',             Maps.HeatMap);
assignin(wsc, 'spMap',             Maps.DispMap);
assignin(wsc, 'spMap_cnt',         Maps.DispCount);
assignin(wsc, 'raw_heatmap',       Maps.RawHeat);
assignin(wsc, 'raw_spatialDisp',   Maps.RawDisp);
assignin(wsc, 'cell_outline',      Maps.CellOutline);
assignin(wsc, 'nuc_outline',       Maps.NucOutlineHeat);
assignin(wsc, 'nuc_contour',       Maps.NucOutlineDisp);
assignin(wsc, 'cntr_low',          Maps.NucCntrBounds(1));
assignin(wsc, 'cntr_up',           Maps.NucCntrBounds(2));
assignin(wsc, 'Simu_spMap',        Maps.SimDispMap);
assignin(wsc, 'Simu_spatialDisp',  Maps.SimSpatDisp);
clear Maps

% -------------------------------------------------------------------------
% 5. Geometry reconstruction ----------------------------------------------
assignin(wsc, 'Nuc_ratio_spMap',  Geometry.NucToCellRatio);
assignin(wsc, 'cell_para',        Geometry.CellParams);
assignin(wsc, 'nuc_para',         Geometry.NucParams);
assignin(wsc, 'Ro',               Geometry.R0);
assignin(wsc, 'Lo',               Geometry.L0);
assignin(wsc, 'convert_factor',   Geometry.ConvertFactor);
assignin(wsc, 'scaling_width',    Geometry.WidthScale);
assignin(wsc, 'scaling_length',   Geometry.LengthScale);
assignin(wsc, 'scale_bar',        Geometry.ScaleBar);
assignin(wsc, 'anchor_htMap',     Geometry.AnchorHeat);
assignin(wsc, 'kb_nuc_htMap',     Geometry.NucKB_Heat);
assignin(wsc, 'anchor_spMap',     Geometry.AnchorDisp);
assignin(wsc, 'kb_nuc_spMap',     Geometry.NucKB_Disp);
clear Geometry

% -------------------------------------------------------------------------
% 6. Volume & accessibility ----------------------------------------------
assignin(wsc, 'lb',               VolAccess.mapBounds(1));
assignin(wsc, 'rb',               VolAccess.mapBounds(2));
assignin(wsc, 'acc_range',        VolAccess.accRange);
assignin(wsc, 'cmp_acc_htMap',    VolAccess.SimuHeatmap);
assignin(wsc, 'nuc_volume_ratio', VolAccess.NucVolRatio);
assignin(wsc, 'acc_fit',          VolAccess.Accessibility);
assignin(wsc, 'Metric_optz',      VolAccess.Similarity);
clear VolAccess

% -------------------------------------------------------------------------
% 7. Brownian‑dynamics simulation -----------------------------------------
assignin(wsc, 'Num_run',          BDSim.NumRuns);
assignin(wsc, 'StepSize_grid',    BDSim.Step);
assignin(wsc, 'Gap_grid',         BDSim.Gap);
assignin(wsc, 'score_w1',         BDSim.WeightScore(1));
assignin(wsc, 'score_w2',         BDSim.WeightScore(2));
assignin(wsc, 'num_agl_seg',      BDSim.angleAugment);
assignin(wsc, 'coord_in_nuc',     BDSim.CoordInNuc);
assignin(wsc, 'coord_in_cyto',    BDSim.CoordInCyto);
assignin(wsc, 'coord_out_cell',   BDSim.CoordOutCell);
assignin(wsc, 'weight_matrix1',   BDSim.WeightMat1);
assignin(wsc, 'weight_matrix2',   BDSim.WeightMat2);
assignin(wsc, 'Store_similarity', BDSim.Similarity);
assignin(wsc, 'Simu_Traj_final',  BDSim.TrajFinal);
assignin(wsc, 'sel_stepsize',     BDSim.SelStepGap(1));
assignin(wsc, 'sel_gap',          BDSim.SelStepGap(2));
clear BDSim

% -------------------------------------------------------------------------
% 8. Raw masks, fits, tracks ---------------------------------------------
if strcmpi(rawData,'on')
    load(fullfile(result_folder,'RawData.mat'),'RawData');
    assignin(wsc, 'Store_mask',  RawData.Mask);
    assignin(wsc, 'Store_fits',  RawData.Fit);
    assignin(wsc, 'Store_track', RawData.Track);
    clear RawData
end
fprintf('\n All variables restored to workspace.\n');
end
=======
function UnpackResults(result_folder, rawData)

% -------------------------------------------------------------------------
% 0. Locate results --------------------------------------------------------
% -------------------------------------------------------------------------
if nargin < 1 || isempty(result_folder)
    result_folder = uigetdir(pwd, 'Select result folder');
    if result_folder == 0
        disp('No folder was selected');
        return
    end
end

if nargin < 2 || isempty(rawData)
    rawData = 'off';
end


% -------------------------------------------------------------------------
% 1. Load structures -------------------------------------------------------
% -------------------------------------------------------------------------
load(fullfile(result_folder, 'Results.mat'), ...
     'LocFit','CellMorph','Maps','Geometry','VolAccess','BDSim');

% Put everything into the caller workspace so the variables behave as if
% they had just been computed there.
wsc = 'caller';

% -------------------------------------------------------------------------
% 2. LocFit ---------------------------------------------------------------
assignin(wsc, 'pixel_size',  LocFit.PixelSize);
assignin(wsc, 'locFitCI',    LocFit.Conf95);
assignin(wsc, 'locFitCI_UB', LocFit.CI_UB);
clear LocFit

% -------------------------------------------------------------------------
% 3. CellMorph & filtering info -------------------------------------------
assignin(wsc, 'mean_ar',         CellMorph.AspectRatio);
assignin(wsc, 'cell_avg_width',  CellMorph.avgWidth);
assignin(wsc, 'morph_by_mov',    CellMorph.MorphByMovie);
assignin(wsc, 'sel_filter',      CellMorph.Filter);
assignin(wsc, 'low_bd',          CellMorph.filterBounds(1));
assignin(wsc, 'up_bd',           CellMorph.filterBounds(2));
assignin(wsc, 'cellIndices',     CellMorph.ValidIdx);
clear CellMorph

% -------------------------------------------------------------------------
% 4. Heat‑ & displacement maps --------------------------------------------
assignin(wsc, 'aspt_ratio',        Maps.CellAR);
assignin(wsc, 'num_pix_y',         Maps.numYGrid);
assignin(wsc, 'htMap',             Maps.HeatMap);
assignin(wsc, 'spMap',             Maps.DispMap);
assignin(wsc, 'spMap_cnt',         Maps.DispCount);
assignin(wsc, 'raw_heatmap',       Maps.RawHeat);
assignin(wsc, 'raw_spatialDisp',   Maps.RawDisp);
assignin(wsc, 'cell_outline',      Maps.CellOutline);
assignin(wsc, 'nuc_outline',       Maps.NucOutlineHeat);
assignin(wsc, 'nuc_contour',       Maps.NucOutlineDisp);
assignin(wsc, 'cntr_low',          Maps.NucCntrBounds(1));
assignin(wsc, 'cntr_up',           Maps.NucCntrBounds(2));
assignin(wsc, 'Simu_spMap',        Maps.SimDispMap);
assignin(wsc, 'Simu_spatialDisp',  Maps.SimSpatDisp);
clear Maps

% -------------------------------------------------------------------------
% 5. Geometry reconstruction ----------------------------------------------
assignin(wsc, 'Nuc_ratio_spMap',  Geometry.NucToCellRatio);
assignin(wsc, 'cell_para',        Geometry.CellParams);
assignin(wsc, 'nuc_para',         Geometry.NucParams);
assignin(wsc, 'Ro',               Geometry.R0);
assignin(wsc, 'Lo',               Geometry.L0);
assignin(wsc, 'convert_factor',   Geometry.ConvertFactor);
assignin(wsc, 'scaling_width',    Geometry.WidthScale);
assignin(wsc, 'scaling_length',   Geometry.LengthScale);
assignin(wsc, 'scale_bar',        Geometry.ScaleBar);
assignin(wsc, 'anchor_htMap',     Geometry.AnchorHeat);
assignin(wsc, 'kb_nuc_htMap',     Geometry.NucKB_Heat);
assignin(wsc, 'anchor_spMap',     Geometry.AnchorDisp);
assignin(wsc, 'kb_nuc_spMap',     Geometry.NucKB_Disp);
clear Geometry

% -------------------------------------------------------------------------
% 6. Volume & accessibility ----------------------------------------------
assignin(wsc, 'lb',               VolAccess.mapBounds(1));
assignin(wsc, 'rb',               VolAccess.mapBounds(2));
assignin(wsc, 'acc_range',        VolAccess.accRange);
assignin(wsc, 'cmp_acc_htMap',    VolAccess.SimuHeatmap);
assignin(wsc, 'nuc_volume_ratio', VolAccess.NucVolRatio);
assignin(wsc, 'acc_fit',          VolAccess.Accessibility);
assignin(wsc, 'Metric_optz',      VolAccess.Similarity);
clear VolAccess

% -------------------------------------------------------------------------
% 7. Brownian‑dynamics simulation -----------------------------------------
assignin(wsc, 'Num_run',          BDSim.NumRuns);
assignin(wsc, 'StepSize_grid',    BDSim.Step);
assignin(wsc, 'Gap_grid',         BDSim.Gap);
assignin(wsc, 'score_w1',         BDSim.WeightScore(1));
assignin(wsc, 'score_w2',         BDSim.WeightScore(2));
assignin(wsc, 'num_agl_seg',      BDSim.angleAugment);
assignin(wsc, 'coord_in_nuc',     BDSim.CoordInNuc);
assignin(wsc, 'coord_in_cyto',    BDSim.CoordInCyto);
assignin(wsc, 'coord_out_cell',   BDSim.CoordOutCell);
assignin(wsc, 'weight_matrix1',   BDSim.WeightMat1);
assignin(wsc, 'weight_matrix2',   BDSim.WeightMat2);
assignin(wsc, 'Store_similarity', BDSim.Similarity);
assignin(wsc, 'Simu_Traj_final',  BDSim.TrajFinal);
assignin(wsc, 'sel_stepsize',     BDSim.SelStepGap(1));
assignin(wsc, 'sel_gap',          BDSim.SelStepGap(2));
clear BDSim

% -------------------------------------------------------------------------
% 8. Raw masks, fits, tracks ---------------------------------------------
if strcmpi(rawData,'on')
    load(fullfile(result_folder,'RawData.mat'),'RawData');
    assignin(wsc, 'Store_mask',  RawData.Mask);
    assignin(wsc, 'Store_fits',  RawData.Fit);
    assignin(wsc, 'Store_track', RawData.Track);
    clear RawData
end
fprintf('\n All variables restored to workspace.\n');
end
>>>>>>> 5e4cf7b (Initial commit)
