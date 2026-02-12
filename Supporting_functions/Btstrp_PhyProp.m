<<<<<<< HEAD
function [Btstrp_NucVol, Btstrp_Acc, Btstrp_Vis_Nuc, ...
          Btstrp_Vis_Cyto, Btstrp_Vis_PeriCore, Btstrp_Vis_Domain,...
          Best_hyperpara] = ...
          Btstrp_PhyProp( ...
          bootstrp_list,   Store_track,  Store_mask,       Store_fits,  ...
          morph_by_mov,    num_pix_y,    locPrecisionUB,   pixel_size,     cell_avg_width, ...
          cmp_range,       acc_metric,   Num_simu,         StepSize_grid,  Gap_grid, ...
          sp_weight,       cntr_margin,  MacroDomain,      score_weight,   constant_setup, ...
          slack_Nuc,       slack_Peri,   slack_Core,       slack_Domain)
%--------------------------------------------------------------------------
%  Btstrp_PhyProp : Bootstrapping with cells for calculations of nucleoid volume, 
%                   accessibility and viscosities in rod-shaped bacterial cells.
%
%  INPUTS (abbreviated list)
%    bootstrp_list : cell array where each cell holds indices for one
%                    bootstrap replicate.
%    Store_*       : cell arrays with tracking data, masks, fits, ...
%    morph_by_mov  : either a scalar aspect ratio or a cell array that
%                    supplies per-movie cell morphology information.
%    num_pix_y     : pixel count in the short axis of the cell patch.
%   locPrecisionUB : Localizaiton preicision upper bound.
%    pixel_size    : pixel size from camera detector.
%    cell_avg_width: average cell width for scalbar calculation (3D cell
%                    gerometry model).
%    cmp_range     : 'all' (default) or 'adjust', specifies region to be
%                    compared when computing heat-map similarity.
%    acc_metric    : Metric used for heatmap comparison to determine
%                    accessibility.
%    Num_simu      : Number of BD simulation for each condition (geometry 
%                    and step size).
%    StepSize_grid : Range of step sizes for BD simulation of each cell
%                    geometry.
%    Gap_grid      : Range of gaps from simulated trajectories to plot
%                    spMap from simulation for similarity evaluation.
%    sp_weight     : Weight scale for spMap weight matrix; (1st 
%                    transformation for L2 Norm, 2nd transformation for PCC).
%                    organized in cell.{[0,1,0],[1,1,0]}
%    cntr_margin   : Contour percentage margin to determin boubdary between
%                    nucleoid periphery and core.
%    MacroDomain   : Macrodomain outlines determined by Locus tracking
%                    experiments. Organized in cell for further usage.
%    score_weight  : weight for L2 Norm score and PCC store in the spMap
%                    comparison. default [0.5, -100].
%   constant_setup : [delta_t,temp_K,probe_radius] in ms, K, and m.
%    slack_*       : range for parameters used in 2-pop Rayleigh distribution 
%                    fitting. [slack_d, slack_acc].
%
%  OUTPUTS (NaN if unavailable)
%    Btstrp_NucVol        : bootstrap result of nucleoid volume ratio
%    Btstrp_Acc           : bootstrap result of nucleoid accessibility
%    Btstrp_Vis_Nuc       : bootstrap result of nucleoidf viscosity 
%    Btstrp_Vis_Cyto      : bootstrap result of cytoplasm viscosity 
%    Btstrp_Vis_PeriCore  : [Vis_periphery  Vis_core]  
%    Btstrp_Vis_Domain    : viscosity for each macrodomain 
%    Best_hyperpara       : optimal step size and gap for simulated spMap
%
% Author: Xiaofeng Dai
% Date: 06/30/2025
%--------------------------------------------------------------------------

%% ------------------------------------------------------------------------
%% 0. Defensive argument handling
%-------------------------------------------------------------------------%
% bootstrp_list must be cell array
if ~iscell(bootstrp_list)
    error('bootstrp_list must be a cell array.');
end

% The three “Store_*” cell arrays must have identical sizes
if ~(isequal(size(Store_track), size(Store_mask)) && ...
     isequal(size(Store_track), size(Store_fits)))
    error('Store_track, Store_mask and Store_fits must be the same size.');
end

% morph_by_mov: either a scalar or a cell array with the same size
if ~(isscalar(morph_by_mov) && isnumeric(morph_by_mov)) && ...
   ~(iscell(morph_by_mov)  && isequal(size(morph_by_mov), size(Store_track)))
    error(['morph_by_mov must be either a numeric scalar or a cell array ' ...
           'with the same size as Store_track/Store_mask/Store_fits.']);
end


% cmp_range can only be 'all' or 'adjust'
if ~strcmpi(cmp_range,'all') && ~strcmpi(cmp_range,'adjust')
    warning('cmp_range "%s" not recognised; defaulting to "all".',cmp_range);
    cmp_range = 'all';
end

% acc_metric has to be chosen from L2 Norm, PCC similarity, SSIM or
% Combine
if ~ismember(acc_metric, {'L2 Norm','PCC','SSIM','Combine'})
    error('wrong metric for heatmap comparison.');
end

if nargin < 18 || isempty(score_weight)
    score_weight = [0.5, -100]; 
end

% constant_setup default as 40ms, 303K, 12.5 nm as probe radius
if nargin < 19 || isempty(constant_setup)
    constant_setup = [40, 303, 12.5e-9];     
elseif numel(constant_setup) ~= 3
    error('constant_setup must be a 1×3 numeric vector.');
end

% default setting for slack_*
if nargin < 20 || isempty(slack_Nuc)
    slack_Nuc = [0.15, 0.075];                
elseif numel(slack_Nuc) ~= 2
    error('slack_Nuc must be a 1×2 numeric vector.');
end

if nargin < 21 || isempty(slack_Peri)
    slack_Peri = [0.2, 0.1];                
elseif numel(slack_Peri) ~= 2
    error('slack_Peri must be a 1×2 numeric vector.');
end

if nargin < 22 || isempty(slack_Core)
    slack_Core = [0.2, 0.1];              
elseif numel(slack_Core) ~= 2
    error('slack_Core must be a 1×2 numeric vector.');
end

if nargin < 23 || isempty(slack_Domain)
    slack_Domain = [0.2, 0.1];             
elseif numel(slack_Domain) ~= 2
    error('slack_Domain must be a 1×2 numeric vector.');
end

nBoot = numel(bootstrp_list);
acc_range = linspace(0,1,101);                   

%% ------------------------------------------------------------------------
%% 1. Pre-allocate outputs
%-------------------------------------------------------------------------%
Btstrp_NucVol         = NaN(nBoot,1);
Btstrp_Acc            = NaN(nBoot,1);
Btstrp_Vis_Nuc        = NaN(nBoot,1);
Btstrp_Vis_Cyto       = NaN(nBoot,1);
Btstrp_Vis_PeriCore   = NaN(nBoot,2);           % [peri  core]
Best_hyperpara        = NaN(nBoot,2);           % [stepsize, gap]

if ~isempty(MacroDomain)
    nDom                   = numel(MacroDomain);
    Btstrp_Vis_Domain      = NaN(nBoot,nDom);
else
    Btstrp_Vis_Domain      = [];              
end

%% ------------------------------------------------------------------------
%% 2. Loop over bootstrap replicates  (PARFOR for speed)
%-------------------------------------------------------------------------%
locPrecPix = locPrecisionUB / pixel_size;        % localisation precision in px

parfor ii = 1:nBoot        
    msg = '';                                    % diagnostic text
    try
        %% 2.1 Re-organise index list
        idxCell = reorgIndex( bootstrp_list{ii} );  

        %% 2.2 Determine the aspect ratio (long/short axis)
        aspt_ratio = determineAspectRatio( ...
                        morph_by_mov, idxCell, num_pix_y);

        %% 2.3 Build heat map and spatial displacement map
        [raw_ht, htMap] = Plot_heatmap_withFilter( ...
            Store_track, Store_mask, Store_fits, idxCell, ...
            1/num_pix_y, aspt_ratio, locPrecPix, ...
            'off','on',slanCM('viridis'));          

        [raw_disp, spMap, ~, spMap_cnt] = Spatial_Displacement_withFilter( ...
            Store_track, Store_mask, Store_fits, idxCell, pixel_size, ...
            3, 1/num_pix_y, aspt_ratio, 0.05, ...
            locPrecPix,'off','on',slanCM('torch'));

        %% 2.4 Cell & nucleoid outlines
        try
            [cell_outline,nuc_outline] = Outline_by_htMap( ...
                htMap, raw_ht, 'interpolation','gradient', ...
                aspt_ratio, num_pix_y, 'off','on', 10,'off');
        catch
            msg = sprintf('%s|Using KDE fallback in replicate %d',msg,ii);
            [cell_outline,nuc_outline] = Outline_by_htMap( ...
                htMap, raw_ht, 'kde','gradient', ...
                aspt_ratio, num_pix_y, 'off','on', 5,'off');
        end

        %% 2.5 Outline from spatial map – returns two contours
        [nuc_contours, ~, cntr_low] = ...
                Outline_by_spMap(spMap, nuc_outline,'off',cell_outline, ...
                                  aspt_ratio, num_pix_y);

        %% 2.6 Geometric fit of the cell outline
        [Rfit,Lfit] = Rod_Morph_fit(cell_outline,num_pix_y,aspt_ratio,'off');

        %---   rescale to a canonical rod  (radius = 5 arbitrary units)   --%
        R0 = 5;
        cFactor        = R0 / Rfit;              % pixel->model scale
        L0             = round(Lfit * cFactor * 10) / 10;
        scalePixShort  = num_pix_y * cFactor / 2;
        scalePixLong   = num_pix_y * aspt_ratio * cFactor / 2;
        scaleBar_nm    = cell_avg_width/(2*scalePixShort)*1000; % 1 model-unit in nm

        %% 2.7 Decide inner vs outer nucleoid contour
        useInner =   max(nuc_contours{1})./[num_pix_y*aspt_ratio, num_pix_y] > [0.915 0.85]; %[0.875 0.8]
        if all(useInner)
            msg = sprintf('%s|Outer boundary dubious – Use inner boundary.',msg);
            selContour  = nuc_contours{2};
        else
            selContour  = nuc_contours{1};
        end
        [anchor_sp, kb_nuc_sp] = NucAnchor( ...
            selContour, aspt_ratio, cFactor, num_pix_y,'off');

        %% 2.8 Choose comparison window in x (long axis)
        [lb_x, rb_x] = comparisonWindow(cmp_range, nuc_outline, ...
                                        num_pix_y, aspt_ratio);

        %% 2.9 STATIC simulation – accessibility
        [nucVolRatio, cmpAccStack] = Nuc_acc( ...
            2e6, acc_range, R0, L0, anchor_sp, ...
            scalePixShort, scalePixLong, num_pix_y, aspt_ratio, ...
            2*size(raw_ht,1), 'off');

        Metric = htMap_similarity( ...
                    htMap, cmpAccStack, acc_range, lb_x, rb_x,'off');

        % Get accessibility from similarity quantification
        acc_fit = AccFit0(acc_range,Metric,acc_metric);
        

        %% 2.10 DYNAMIC simulation – single-particle tracks replay
        [simDisp,coord_in_nuc,coord_in_cyto, ~, ~, sel_stepsize, sel_gap, ~] = DynamicSimulation( ...
            Num_simu, StepSize_grid, Gap_grid, spMap, spMap_cnt, sp_weight, ...
            cell_outline, anchor_sp, num_pix_y, aspt_ratio, cFactor, ...
            scaleBar_nm, scalePixShort, scalePixLong, R0, L0, ...
            kb_nuc_sp, score_weight);

        %% 2.11 Physical constants
        delta_t       = constant_setup(1);    % ms
        temp_K        = constant_setup(2);    % K
        probe_radius  = constant_setup(3);    % nm

        %% 2.12 Viscosity in nucleoid vs cytosol
        [visNuc, ~, ~, visCyt, ~, ~] = CalVis_Nuc_Cyto( ...
            delta_t, temp_K, probe_radius, acc_fit, R0, L0, ...
            anchor_sp, spMap, scalePixShort, scalePixLong, ...
            num_pix_y, aspt_ratio, selContour, raw_disp, simDisp, ...
            coord_in_nuc,coord_in_cyto, slack_Nuc);                      

        %% 2.13  Viscosities in periphery and core  (optional)
        visPeriCore = [NaN NaN];                 
        
        if ~isempty(cntr_margin)
            try
                [vPeri,~,~,vCore,~,~] = CalVis_Peri_Core( ...
                    delta_t, temp_K, probe_radius, spMap, ...
                    [cntr_low+cntr_margin, cntr_low], ...
                    cell_outline, num_pix_y, aspt_ratio, ...
                    raw_disp, simDisp, acc_fit, ...
                    R0, L0, anchor_sp, ...
                    scalePixShort, scalePixLong, ...
                    slack_Peri, slack_Core);
        
                visPeriCore = [vPeri, vCore];   
        
            catch ME
                msg = sprintf('%s|Peri/Core Vis fail: %s', msg, ME.message);
                % visPeriCore remains [NaN NaN]
            end
        end

        %% 2.14  Macro-domain viscosities  (optional)
        nDom   = max(1, numel(MacroDomain));
        visDom = NaN(1, nDom);                  
        
        if ~isempty(MacroDomain)
            for kk = 1:numel(MacroDomain)
                try
                    DomOutline = MacroDomain{kk};
                    [vTmp,~,~] = CalVis_MacroDomain( ...
                        delta_t, temp_K, probe_radius, spMap, DomOutline, ...
                        nuc_outline, num_pix_y, aspt_ratio, ...
                        raw_disp, simDisp, acc_fit, ...
                        R0, L0, anchor_sp, scalePixShort, scalePixLong, ...
                        slack_Domain);
        
                    visDom(kk) = vTmp;          % overwrite only on success
        
                catch ME
                    msg = sprintf('%s|Macrodomain Vis %d fail: %s', msg, kk, ME.message);
                    % leave NaN
                end
            end
        end

        %% 2.15 Save results
        Best_hyperpara(ii,:)    = [sel_stepsize, sel_gap];
        Btstrp_NucVol(ii)       = nucVolRatio;
        Btstrp_Acc(ii)          = acc_fit;
        Btstrp_Vis_Nuc(ii)      = visNuc;
        Btstrp_Vis_Cyto(ii)     = visCyt;
        Btstrp_Vis_PeriCore(ii,:) = visPeriCore;
        if ~isempty(MacroDomain);  Btstrp_Vis_Domain(ii,:) = visDom;  end

    catch ME
        % record first frame (where the error really happened)
        if ~isempty(ME.stack)
            stk = ME.stack(1);
            msg = sprintf('%s|%s (%s:%d)', ...
                  msg, ME.message, stk.name, stk.line);
        else
            msg = sprintf('%s|%s', msg, ME.message);
        end
    end

    % Print one concise line per replicate (outside try so it never fails)
    if ~isempty(msg);  fprintf('[%3d/%3d] %s\n',ii,nBoot,msg);  end
end   % parfor

end  % Btstrp_PhyProp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Local helper functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function aspt_ratio = determineAspectRatio(morph_by_mov, idxCell, num_pix_y)
% Return aspect ratio for the current bootstrap replicate
    if isscalar(morph_by_mov) && isnumeric(morph_by_mov)
        aspt_ratio = morph_by_mov;

    elseif iscell(morph_by_mov) && numel(morph_by_mov)==numel(idxCell)
        % Iterate through all movies in the replicate and collect the
        % aspect ratios of the selected cells
        aspVals = cell(numel(idxCell),1);
        for ia = 1:numel(idxCell)
            subIdx   = idxCell{ia};
            morphTbl = morph_by_mov{ia};      % user supplied table/array
            asp      = nan(numel(subIdx),1);
            for ib = 1:numel(subIdx)
                id       = subIdx(ib);
                rowMask  = morphTbl(:,1)==id;
                asp(ib)  = morphTbl(rowMask,4);
            end
            aspVals{ia} = asp;
        end
        aspt_ratio = round( mean(cat(1,aspVals{:})) * num_pix_y ) / num_pix_y;

    else
        error('morph_by_mov has incompatible type/size.');
    end
end
%-------------------------------------------------------------------------%
function [lb, rb] = comparisonWindow(mode, nuc_outline, nPixY, aspt_ratio)
% Pick x-range used for heat-map similarity
    switch lower(mode)
        case 'all'
            lb = 1;
            rb = nPixY * aspt_ratio;
        case 'adjust'
            nPts   = size(nuc_outline,1);
            nPick  = max(1, round(nPts*0.05));
            rb = floor(mean(maxk(nuc_outline(:,1),nPick)));
            lb = floor(mean(mink(nuc_outline(:,1),nPick)));
        otherwise
            lb = 1;
            rb = nPixY * aspt_ratio;
    end
end
%-------------------------------------------------------------------------%
function acc_fit = AccFit0(acc_range,Metric_optz,method)
% Calculate accessibility based on similarity scores from Metric_optz
    acc_range_fine = linspace(0,1,501);
    acc_fit_all = NaN(size(Metric_optz,2),1);
    for i = 1:size(Metric_optz,2)
        f = fit(acc_range.',Metric_optz(:,i),'poly5');
        if i ==1
            [~,acc_idx] = min(f(linspace(0,1,501)));
        else
            [~,acc_idx] = max(f(linspace(0,1,501)));
        end
        acc_fit_all(i) = acc_range_fine(acc_idx);
    end
    switch method
        case 'L2 Norm'
            acc_fit = acc_fit_all(1);
        case 'PCC'
            acc_fit = acc_fit_all(2);
        case 'SSIM'
            acc_fit = acc_fit_all(3);
        case 'Combine'
            acc_fit = mean(rmoutliers(acc_fit_all));
        otherwise
            disp('Wrong method, Combine method is used')
            acc_fit = mean(rmoutliers(acc_fit_all));
    end
    % fprintf('Nucleoid Accessibility = %.4f\n', acc_fit);
end
=======
function [Btstrp_NucVol, Btstrp_Acc, Btstrp_Vis_Nuc, ...
          Btstrp_Vis_Cyto, Btstrp_Vis_PeriCore, Btstrp_Vis_Domain,...
          Best_hyperpara] = ...
          Btstrp_PhyProp( ...
          bootstrp_list,   Store_track,  Store_mask,       Store_fits,  ...
          morph_by_mov,    num_pix_y,    locPrecisionUB,   pixel_size,     cell_avg_width, ...
          cmp_range,       acc_metric,   Num_simu,         StepSize_grid,  Gap_grid, ...
          sp_weight,       cntr_margin,  MacroDomain,      score_weight,   constant_setup, ...
          slack_Nuc,       slack_Peri,   slack_Core,       slack_Domain)
%--------------------------------------------------------------------------
%  Btstrp_PhyProp : Bootstrapping with cells for calculations of nucleoid volume, 
%                   accessibility and viscosities in rod-shaped bacterial cells.
%
%  INPUTS (abbreviated list)
%    bootstrp_list : cell array where each cell holds indices for one
%                    bootstrap replicate.
%    Store_*       : cell arrays with tracking data, masks, fits, ...
%    morph_by_mov  : either a scalar aspect ratio or a cell array that
%                    supplies per-movie cell morphology information.
%    num_pix_y     : pixel count in the short axis of the cell patch.
%   locPrecisionUB : Localizaiton preicision upper bound.
%    pixel_size    : pixel size from camera detector.
%    cell_avg_width: average cell width for scalbar calculation (3D cell
%                    gerometry model).
%    cmp_range     : 'all' (default) or 'adjust', specifies region to be
%                    compared when computing heat-map similarity.
%    acc_metric    : Metric used for heatmap comparison to determine
%                    accessibility.
%    Num_simu      : Number of BD simulation for each condition (geometry 
%                    and step size).
%    StepSize_grid : Range of step sizes for BD simulation of each cell
%                    geometry.
%    Gap_grid      : Range of gaps from simulated trajectories to plot
%                    spMap from simulation for similarity evaluation.
%    sp_weight     : Weight scale for spMap weight matrix; (1st 
%                    transformation for L2 Norm, 2nd transformation for PCC).
%                    organized in cell.{[0,1,0],[1,1,0]}
%    cntr_margin   : Contour percentage margin to determin boubdary between
%                    nucleoid periphery and core.
%    MacroDomain   : Macrodomain outlines determined by Locus tracking
%                    experiments. Organized in cell for further usage.
%    score_weight  : weight for L2 Norm score and PCC store in the spMap
%                    comparison. default [0.5, -100].
%   constant_setup : [delta_t,temp_K,probe_radius] in ms, K, and m.
%    slack_*       : range for parameters used in 2-pop Rayleigh distribution 
%                    fitting. [slack_d, slack_acc].
%
%  OUTPUTS (NaN if unavailable)
%    Btstrp_NucVol        : bootstrap result of nucleoid volume ratio
%    Btstrp_Acc           : bootstrap result of nucleoid accessibility
%    Btstrp_Vis_Nuc       : bootstrap result of nucleoidf viscosity 
%    Btstrp_Vis_Cyto      : bootstrap result of cytoplasm viscosity 
%    Btstrp_Vis_PeriCore  : [Vis_periphery  Vis_core]  
%    Btstrp_Vis_Domain    : viscosity for each macrodomain 
%    Best_hyperpara       : optimal step size and gap for simulated spMap
%
% Author: Xiaofeng Dai
% Date: 06/30/2025
%--------------------------------------------------------------------------

%% ------------------------------------------------------------------------
%% 0. Defensive argument handling
%-------------------------------------------------------------------------%
% bootstrp_list must be cell array
if ~iscell(bootstrp_list)
    error('bootstrp_list must be a cell array.');
end

% The three “Store_*” cell arrays must have identical sizes
if ~(isequal(size(Store_track), size(Store_mask)) && ...
     isequal(size(Store_track), size(Store_fits)))
    error('Store_track, Store_mask and Store_fits must be the same size.');
end

% morph_by_mov: either a scalar or a cell array with the same size
if ~(isscalar(morph_by_mov) && isnumeric(morph_by_mov)) && ...
   ~(iscell(morph_by_mov)  && isequal(size(morph_by_mov), size(Store_track)))
    error(['morph_by_mov must be either a numeric scalar or a cell array ' ...
           'with the same size as Store_track/Store_mask/Store_fits.']);
end


% cmp_range can only be 'all' or 'adjust'
if ~strcmpi(cmp_range,'all') && ~strcmpi(cmp_range,'adjust')
    warning('cmp_range "%s" not recognised; defaulting to "all".',cmp_range);
    cmp_range = 'all';
end

% acc_metric has to be chosen from L2 Norm, PCC similarity, SSIM or
% Combine
if ~ismember(acc_metric, {'L2 Norm','PCC','SSIM','Combine'})
    error('wrong metric for heatmap comparison.');
end

if nargin < 18 || isempty(score_weight)
    score_weight = [0.5, -100]; 
end

% constant_setup default as 40ms, 303K, 12.5 nm as probe radius
if nargin < 19 || isempty(constant_setup)
    constant_setup = [40, 303, 12.5e-9];     
elseif numel(constant_setup) ~= 3
    error('constant_setup must be a 1×3 numeric vector.');
end

% default setting for slack_*
if nargin < 20 || isempty(slack_Nuc)
    slack_Nuc = [0.15, 0.075];                
elseif numel(slack_Nuc) ~= 2
    error('slack_Nuc must be a 1×2 numeric vector.');
end

if nargin < 21 || isempty(slack_Peri)
    slack_Peri = [0.2, 0.1];                
elseif numel(slack_Peri) ~= 2
    error('slack_Peri must be a 1×2 numeric vector.');
end

if nargin < 22 || isempty(slack_Core)
    slack_Core = [0.2, 0.1];              
elseif numel(slack_Core) ~= 2
    error('slack_Core must be a 1×2 numeric vector.');
end

if nargin < 23 || isempty(slack_Domain)
    slack_Domain = [0.2, 0.1];             
elseif numel(slack_Domain) ~= 2
    error('slack_Domain must be a 1×2 numeric vector.');
end

nBoot = numel(bootstrp_list);
acc_range = linspace(0,1,101);                   

%% ------------------------------------------------------------------------
%% 1. Pre-allocate outputs
%-------------------------------------------------------------------------%
Btstrp_NucVol         = NaN(nBoot,1);
Btstrp_Acc            = NaN(nBoot,1);
Btstrp_Vis_Nuc        = NaN(nBoot,1);
Btstrp_Vis_Cyto       = NaN(nBoot,1);
Btstrp_Vis_PeriCore   = NaN(nBoot,2);           % [peri  core]
Best_hyperpara        = NaN(nBoot,2);           % [stepsize, gap]

if ~isempty(MacroDomain)
    nDom                   = numel(MacroDomain);
    Btstrp_Vis_Domain      = NaN(nBoot,nDom);
else
    Btstrp_Vis_Domain      = [];              
end

%% ------------------------------------------------------------------------
%% 2. Loop over bootstrap replicates  (PARFOR for speed)
%-------------------------------------------------------------------------%
locPrecPix = locPrecisionUB / pixel_size;        % localisation precision in px

parfor ii = 1:nBoot        
    msg = '';                                    % diagnostic text
    try
        %% 2.1 Re-organise index list
        idxCell = reorgIndex( bootstrp_list{ii} );  

        %% 2.2 Determine the aspect ratio (long/short axis)
        aspt_ratio = determineAspectRatio( ...
                        morph_by_mov, idxCell, num_pix_y);

        %% 2.3 Build heat map and spatial displacement map
        [raw_ht, htMap] = Plot_heatmap_withFilter( ...
            Store_track, Store_mask, Store_fits, idxCell, ...
            1/num_pix_y, aspt_ratio, locPrecPix, ...
            'off','on',slanCM('viridis'));          

        [raw_disp, spMap, ~, spMap_cnt] = Spatial_Displacement_withFilter( ...
            Store_track, Store_mask, Store_fits, idxCell, pixel_size, ...
            3, 1/num_pix_y, aspt_ratio, 0.05, ...
            locPrecPix,'off','on',slanCM('torch'));

        %% 2.4 Cell & nucleoid outlines
        try
            [cell_outline,nuc_outline] = Outline_by_htMap( ...
                htMap, raw_ht, 'interpolation','gradient', ...
                aspt_ratio, num_pix_y, 'off','on', 10,'off');
        catch
            msg = sprintf('%s|Using KDE fallback in replicate %d',msg,ii);
            [cell_outline,nuc_outline] = Outline_by_htMap( ...
                htMap, raw_ht, 'kde','gradient', ...
                aspt_ratio, num_pix_y, 'off','on', 5,'off');
        end

        %% 2.5 Outline from spatial map – returns two contours
        [nuc_contours, ~, cntr_low] = ...
                Outline_by_spMap(spMap, nuc_outline,'off',cell_outline, ...
                                  aspt_ratio, num_pix_y);

        %% 2.6 Geometric fit of the cell outline
        [Rfit,Lfit] = Rod_Morph_fit(cell_outline,num_pix_y,aspt_ratio,'off');

        %---   rescale to a canonical rod  (radius = 5 arbitrary units)   --%
        R0 = 5;
        cFactor        = R0 / Rfit;              % pixel->model scale
        L0             = round(Lfit * cFactor * 10) / 10;
        scalePixShort  = num_pix_y * cFactor / 2;
        scalePixLong   = num_pix_y * aspt_ratio * cFactor / 2;
        scaleBar_nm    = cell_avg_width/(2*scalePixShort)*1000; % 1 model-unit in nm

        %% 2.7 Decide inner vs outer nucleoid contour
        useInner =   max(nuc_contours{1})./[num_pix_y*aspt_ratio, num_pix_y] > [0.915 0.85]; %[0.875 0.8]
        if all(useInner)
            msg = sprintf('%s|Outer boundary dubious – Use inner boundary.',msg);
            selContour  = nuc_contours{2};
        else
            selContour  = nuc_contours{1};
        end
        [anchor_sp, kb_nuc_sp] = NucAnchor( ...
            selContour, aspt_ratio, cFactor, num_pix_y,'off');

        %% 2.8 Choose comparison window in x (long axis)
        [lb_x, rb_x] = comparisonWindow(cmp_range, nuc_outline, ...
                                        num_pix_y, aspt_ratio);

        %% 2.9 STATIC simulation – accessibility
        [nucVolRatio, cmpAccStack] = Nuc_acc( ...
            2e6, acc_range, R0, L0, anchor_sp, ...
            scalePixShort, scalePixLong, num_pix_y, aspt_ratio, ...
            2*size(raw_ht,1), 'off');

        Metric = htMap_similarity( ...
                    htMap, cmpAccStack, acc_range, lb_x, rb_x,'off');

        % Get accessibility from similarity quantification
        acc_fit = AccFit0(acc_range,Metric,acc_metric);
        

        %% 2.10 DYNAMIC simulation – single-particle tracks replay
        [simDisp,coord_in_nuc,coord_in_cyto, ~, ~, sel_stepsize, sel_gap, ~] = DynamicSimulation( ...
            Num_simu, StepSize_grid, Gap_grid, spMap, spMap_cnt, sp_weight, ...
            cell_outline, anchor_sp, num_pix_y, aspt_ratio, cFactor, ...
            scaleBar_nm, scalePixShort, scalePixLong, R0, L0, ...
            kb_nuc_sp, score_weight);

        %% 2.11 Physical constants
        delta_t       = constant_setup(1);    % ms
        temp_K        = constant_setup(2);    % K
        probe_radius  = constant_setup(3);    % nm

        %% 2.12 Viscosity in nucleoid vs cytosol
        [visNuc, ~, ~, visCyt, ~, ~] = CalVis_Nuc_Cyto( ...
            delta_t, temp_K, probe_radius, acc_fit, R0, L0, ...
            anchor_sp, spMap, scalePixShort, scalePixLong, ...
            num_pix_y, aspt_ratio, selContour, raw_disp, simDisp, ...
            coord_in_nuc,coord_in_cyto, slack_Nuc);                      

        %% 2.13  Viscosities in periphery and core  (optional)
        visPeriCore = [NaN NaN];                 
        
        if ~isempty(cntr_margin)
            try
                [vPeri,~,~,vCore,~,~] = CalVis_Peri_Core( ...
                    delta_t, temp_K, probe_radius, spMap, ...
                    [cntr_low+cntr_margin, cntr_low], ...
                    cell_outline, num_pix_y, aspt_ratio, ...
                    raw_disp, simDisp, acc_fit, ...
                    R0, L0, anchor_sp, ...
                    scalePixShort, scalePixLong, ...
                    slack_Peri, slack_Core);
        
                visPeriCore = [vPeri, vCore];   
        
            catch ME
                msg = sprintf('%s|Peri/Core Vis fail: %s', msg, ME.message);
                % visPeriCore remains [NaN NaN]
            end
        end

        %% 2.14  Macro-domain viscosities  (optional)
        nDom   = max(1, numel(MacroDomain));
        visDom = NaN(1, nDom);                  
        
        if ~isempty(MacroDomain)
            for kk = 1:numel(MacroDomain)
                try
                    DomOutline = MacroDomain{kk};
                    [vTmp,~,~] = CalVis_MacroDomain( ...
                        delta_t, temp_K, probe_radius, spMap, DomOutline, ...
                        nuc_outline, num_pix_y, aspt_ratio, ...
                        raw_disp, simDisp, acc_fit, ...
                        R0, L0, anchor_sp, scalePixShort, scalePixLong, ...
                        slack_Domain);
        
                    visDom(kk) = vTmp;          % overwrite only on success
        
                catch ME
                    msg = sprintf('%s|Macrodomain Vis %d fail: %s', msg, kk, ME.message);
                    % leave NaN
                end
            end
        end

        %% 2.15 Save results
        Best_hyperpara(ii,:)    = [sel_stepsize, sel_gap];
        Btstrp_NucVol(ii)       = nucVolRatio;
        Btstrp_Acc(ii)          = acc_fit;
        Btstrp_Vis_Nuc(ii)      = visNuc;
        Btstrp_Vis_Cyto(ii)     = visCyt;
        Btstrp_Vis_PeriCore(ii,:) = visPeriCore;
        if ~isempty(MacroDomain);  Btstrp_Vis_Domain(ii,:) = visDom;  end

    catch ME
        % record first frame (where the error really happened)
        if ~isempty(ME.stack)
            stk = ME.stack(1);
            msg = sprintf('%s|%s (%s:%d)', ...
                  msg, ME.message, stk.name, stk.line);
        else
            msg = sprintf('%s|%s', msg, ME.message);
        end
    end

    % Print one concise line per replicate (outside try so it never fails)
    if ~isempty(msg);  fprintf('[%3d/%3d] %s\n',ii,nBoot,msg);  end
end   % parfor

end  % Btstrp_PhyProp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Local helper functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function aspt_ratio = determineAspectRatio(morph_by_mov, idxCell, num_pix_y)
% Return aspect ratio for the current bootstrap replicate
    if isscalar(morph_by_mov) && isnumeric(morph_by_mov)
        aspt_ratio = morph_by_mov;

    elseif iscell(morph_by_mov) && numel(morph_by_mov)==numel(idxCell)
        % Iterate through all movies in the replicate and collect the
        % aspect ratios of the selected cells
        aspVals = cell(numel(idxCell),1);
        for ia = 1:numel(idxCell)
            subIdx   = idxCell{ia};
            morphTbl = morph_by_mov{ia};      % user supplied table/array
            asp      = nan(numel(subIdx),1);
            for ib = 1:numel(subIdx)
                id       = subIdx(ib);
                rowMask  = morphTbl(:,1)==id;
                asp(ib)  = morphTbl(rowMask,4);
            end
            aspVals{ia} = asp;
        end
        aspt_ratio = round( mean(cat(1,aspVals{:})) * num_pix_y ) / num_pix_y;

    else
        error('morph_by_mov has incompatible type/size.');
    end
end
%-------------------------------------------------------------------------%
function [lb, rb] = comparisonWindow(mode, nuc_outline, nPixY, aspt_ratio)
% Pick x-range used for heat-map similarity
    switch lower(mode)
        case 'all'
            lb = 1;
            rb = nPixY * aspt_ratio;
        case 'adjust'
            nPts   = size(nuc_outline,1);
            nPick  = max(1, round(nPts*0.05));
            rb = floor(mean(maxk(nuc_outline(:,1),nPick)));
            lb = floor(mean(mink(nuc_outline(:,1),nPick)));
        otherwise
            lb = 1;
            rb = nPixY * aspt_ratio;
    end
end
%-------------------------------------------------------------------------%
function acc_fit = AccFit0(acc_range,Metric_optz,method)
% Calculate accessibility based on similarity scores from Metric_optz
    acc_range_fine = linspace(0,1,501);
    acc_fit_all = NaN(size(Metric_optz,2),1);
    for i = 1:size(Metric_optz,2)
        f = fit(acc_range.',Metric_optz(:,i),'poly5');
        if i ==1
            [~,acc_idx] = min(f(linspace(0,1,501)));
        else
            [~,acc_idx] = max(f(linspace(0,1,501)));
        end
        acc_fit_all(i) = acc_range_fine(acc_idx);
    end
    switch method
        case 'L2 Norm'
            acc_fit = acc_fit_all(1);
        case 'PCC'
            acc_fit = acc_fit_all(2);
        case 'SSIM'
            acc_fit = acc_fit_all(3);
        case 'Combine'
            acc_fit = mean(rmoutliers(acc_fit_all));
        otherwise
            disp('Wrong method, Combine method is used')
            acc_fit = mean(rmoutliers(acc_fit_all));
    end
    % fprintf('Nucleoid Accessibility = %.4f\n', acc_fit);
end
>>>>>>> 5e4cf7b (Initial commit)
