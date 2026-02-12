<<<<<<< HEAD
function [Simu_spatialDisp,coord_in_nuc,coord_in_cyto, Simu_spMap, Simu_Traj_final, sel_stepsize, sel_gap, Store_similarity] =...
    DynamicSimulation(...
    Num_run,StepSize_grid,Gap_grid,spMap,spMap_cnt,sp_weight,...
    cell_outline,anchor_spMap,num_pix_y,aspt_ratio,convert_factor,...
    scale_bar,scaling_width,scaling_length,Ro,Lo,kb_nuc_spMap,score_weight)

% segment map into nucloid, cytoplasm and outside of cell
[coord_in_nuc,coord_in_cyto,coord_out_cell] =...
    Map_Segment(spMap,cell_outline,anchor_spMap,...
    num_pix_y,aspt_ratio,convert_factor,'off');

% weight matrix for spMap comparision
weight_matrix1 = Weight_Cmp_spMap(spMap_cnt,...
    {coord_in_nuc,coord_in_cyto,coord_out_cell},...
    sp_weight{1},[1,1,1],'linear');
% weight matrix for spMap comparision
weight_matrix2 = Weight_Cmp_spMap(spMap_cnt,...
    {coord_in_nuc,coord_in_cyto,coord_out_cell},...
    sp_weight{2},[1,1,1],'linear');


score_w1 = score_weight(1);
score_w2 = score_weight(2);

Store_similarity = cell(length(StepSize_grid),1);
num_agl_seg = 6;
deter_score = [];
Simu_Traj_final = [];

for k = 1:length(StepSize_grid)
    StepSize = StepSize_grid(k);
    Store_result = BD_Simu_batchGen(Num_run,StepSize,[],...
        scale_bar,Ro,Lo,kb_nuc_spMap,anchor_spMap,...
        'on','Collision',...
        floor((max(StepSize_grid)/StepSize)^1.2)*1e5);
    
    buff = cell(1,length(Store_result));
    for j = 1:length(Store_result)
        buff{j} = Store_result{j}(:,1:3);
    end
    Store_result = buff;

    % Data Augmentation
    % rotate coordinates from simulation
    Store_result = SimuTraj_RotAug(Store_result,num_agl_seg);

    Similarity_sp = zeros(length(Gap_grid),6);
    for i = 1:length(Gap_grid)
        % Generate spMap with simulated trajectories 
        % with different Gap for comparison 
        [Simu_spMap,~] = ...
        Spatial_Displacement_Simu(...
        Store_result,Gap_grid(i),num_pix_y,aspt_ratio,...
        scaling_width,scaling_length,scale_bar,0.02,...
        'on','off');
    
        [sp_L2_unweighted, sp_L2, sp_L2_adj,...
        ~, ~, ~] = sp_cmp_similarity(spMap,...
        Simu_spMap, weight_matrix1, coord_in_nuc);
        [~, ~, ~,sp_SSIM, sp_PCC, sp_JSD] =...
            sp_cmp_similarity(spMap, Simu_spMap,...
                              weight_matrix2, coord_in_nuc);
    
        Similarity_sp(i,:) =...
            [sp_L2_unweighted, sp_L2, sp_L2_adj,...
             sp_SSIM, sp_PCC, sp_JSD];
    end
    
    Score = score_w1*Similarity_sp(:,2) +...
            score_w2*Similarity_sp(:,5);

    if isempty(deter_score)
        deter_score = min(Score);
    end
    if deter_score > min(Score)
        deter_score = min(Score);
        Simu_Traj_final = Store_result;
        % min_step = StepSize_grid(k);
        % min_gap = Gap_grid(Score==deter_score);
    end
    Store_similarity{k} = Similarity_sp;
end

similarity_plot = [];
for i = 1:length(StepSize_grid)
    similarity_plot = cat(2,similarity_plot,...
        score_w1*Store_similarity{i}(:,2)+...
        score_w2*Store_similarity{i}(:,5));

    % similarity_plot = cat(2,similarity_plot,...
    %     Store_similarity{i}(:,1));
end
clear i

[sel_gap_idx, sel_step_idx] =...
    find(similarity_plot == min(similarity_plot,[],"all"));

sel_stepsize = StepSize_grid(sel_step_idx);
sel_gap = Gap_grid(sel_gap_idx);

[Simu_spMap,Simu_spatialDisp] =...
    Spatial_Displacement_Simu(Simu_Traj_final,sel_gap,...
    num_pix_y,aspt_ratio,scaling_width,scaling_length,...
    scale_bar,0.05,'on','off');

=======
function [Simu_spatialDisp,coord_in_nuc,coord_in_cyto, Simu_spMap, Simu_Traj_final, sel_stepsize, sel_gap, Store_similarity] =...
    DynamicSimulation(...
    Num_run,StepSize_grid,Gap_grid,spMap,spMap_cnt,sp_weight,...
    cell_outline,anchor_spMap,num_pix_y,aspt_ratio,convert_factor,...
    scale_bar,scaling_width,scaling_length,Ro,Lo,kb_nuc_spMap,score_weight)

% segment map into nucloid, cytoplasm and outside of cell
[coord_in_nuc,coord_in_cyto,coord_out_cell] =...
    Map_Segment(spMap,cell_outline,anchor_spMap,...
    num_pix_y,aspt_ratio,convert_factor,'off');

% weight matrix for spMap comparision
weight_matrix1 = Weight_Cmp_spMap(spMap_cnt,...
    {coord_in_nuc,coord_in_cyto,coord_out_cell},...
    sp_weight{1},[1,1,1],'linear');
% weight matrix for spMap comparision
weight_matrix2 = Weight_Cmp_spMap(spMap_cnt,...
    {coord_in_nuc,coord_in_cyto,coord_out_cell},...
    sp_weight{2},[1,1,1],'linear');


score_w1 = score_weight(1);
score_w2 = score_weight(2);

Store_similarity = cell(length(StepSize_grid),1);
num_agl_seg = 6;
deter_score = [];
Simu_Traj_final = [];

for k = 1:length(StepSize_grid)
    StepSize = StepSize_grid(k);
    Store_result = BD_Simu_batchGen(Num_run,StepSize,[],...
        scale_bar,Ro,Lo,kb_nuc_spMap,anchor_spMap,...
        'on','Collision',...
        floor((max(StepSize_grid)/StepSize)^1.2)*1e5);
    
    buff = cell(1,length(Store_result));
    for j = 1:length(Store_result)
        buff{j} = Store_result{j}(:,1:3);
    end
    Store_result = buff;

    % Data Augmentation
    % rotate coordinates from simulation
    Store_result = SimuTraj_RotAug(Store_result,num_agl_seg);

    Similarity_sp = zeros(length(Gap_grid),6);
    for i = 1:length(Gap_grid)
        % Generate spMap with simulated trajectories 
        % with different Gap for comparison 
        [Simu_spMap,~] = ...
        Spatial_Displacement_Simu(...
        Store_result,Gap_grid(i),num_pix_y,aspt_ratio,...
        scaling_width,scaling_length,scale_bar,0.02,...
        'on','off');
    
        [sp_L2_unweighted, sp_L2, sp_L2_adj,...
        ~, ~, ~] = sp_cmp_similarity(spMap,...
        Simu_spMap, weight_matrix1, coord_in_nuc);
        [~, ~, ~,sp_SSIM, sp_PCC, sp_JSD] =...
            sp_cmp_similarity(spMap, Simu_spMap,...
                              weight_matrix2, coord_in_nuc);
    
        Similarity_sp(i,:) =...
            [sp_L2_unweighted, sp_L2, sp_L2_adj,...
             sp_SSIM, sp_PCC, sp_JSD];
    end
    
    Score = score_w1*Similarity_sp(:,2) +...
            score_w2*Similarity_sp(:,5);

    if isempty(deter_score)
        deter_score = min(Score);
    end
    if deter_score > min(Score)
        deter_score = min(Score);
        Simu_Traj_final = Store_result;
        % min_step = StepSize_grid(k);
        % min_gap = Gap_grid(Score==deter_score);
    end
    Store_similarity{k} = Similarity_sp;
end

similarity_plot = [];
for i = 1:length(StepSize_grid)
    similarity_plot = cat(2,similarity_plot,...
        score_w1*Store_similarity{i}(:,2)+...
        score_w2*Store_similarity{i}(:,5));

    % similarity_plot = cat(2,similarity_plot,...
    %     Store_similarity{i}(:,1));
end
clear i

[sel_gap_idx, sel_step_idx] =...
    find(similarity_plot == min(similarity_plot,[],"all"));

sel_stepsize = StepSize_grid(sel_step_idx);
sel_gap = Gap_grid(sel_gap_idx);

[Simu_spMap,Simu_spatialDisp] =...
    Spatial_Displacement_Simu(Simu_Traj_final,sel_gap,...
    num_pix_y,aspt_ratio,scaling_width,scaling_length,...
    scale_bar,0.05,'on','off');

>>>>>>> 5e4cf7b (Initial commit)
end