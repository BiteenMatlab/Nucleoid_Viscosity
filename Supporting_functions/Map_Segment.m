<<<<<<< HEAD
function [coord_in_nuc,coord_in_cyto,coord_out_cell] = Map_Segment(Map,cell_outline,nuc_outline,num_pix_y,aspt_ratio,convert_factor,fig_plot)

%% set up cell boundary
cell_bndy_coord = cell_outline;
cell_bndy_coord(:,2) = cell_bndy_coord(:,2)-num_pix_y/2-0.5;
cell_bndy_coord(:,1) = cell_bndy_coord(:,1)-num_pix_y*aspt_ratio/2-0.5;
cell_bndy_coord(cell_bndy_coord(:,2)<0,:) = [];
[coord_in_cell, coord_out_cell] = MapCoord_region(Map,cell_bndy_coord,num_pix_y,aspt_ratio,'off');
%% set up nucleoid boundary
nuc_bndy_coord = nuc_outline/convert_factor;
[coord_in_nuc, coord_out_nuc] = MapCoord_region(Map,nuc_bndy_coord,num_pix_y,aspt_ratio,'off');

%% coordinates of cytoplasma pixels
coord_in_cyto = intersect(coord_out_nuc,coord_in_cell,'rows');
% coord_in_nuc --> nucleoid pixels
% coord_out_cell --> pixels outside of cell

%% scatter plot different regions
if strcmp(fig_plot,'on')
    figure
    scatter(coord_in_nuc(:,1),coord_in_nuc(:,2),5,'filled')
    hold on
    scatter(coord_out_cell(:,1),coord_out_cell(:,2),5,'filled')
    hold on
    scatter(coord_in_cyto(:,1),coord_in_cyto(:,2),5,'filled')
    axis equal
    xlim([0 aspt_ratio*num_pix_y])
    ylim([0 num_pix_y])
    % axis off
end
=======
function [coord_in_nuc,coord_in_cyto,coord_out_cell] = Map_Segment(Map,cell_outline,nuc_outline,num_pix_y,aspt_ratio,convert_factor,fig_plot)

%% set up cell boundary
cell_bndy_coord = cell_outline;
cell_bndy_coord(:,2) = cell_bndy_coord(:,2)-num_pix_y/2-0.5;
cell_bndy_coord(:,1) = cell_bndy_coord(:,1)-num_pix_y*aspt_ratio/2-0.5;
cell_bndy_coord(cell_bndy_coord(:,2)<0,:) = [];
[coord_in_cell, coord_out_cell] = MapCoord_region(Map,cell_bndy_coord,num_pix_y,aspt_ratio,'off');
%% set up nucleoid boundary
nuc_bndy_coord = nuc_outline/convert_factor;
[coord_in_nuc, coord_out_nuc] = MapCoord_region(Map,nuc_bndy_coord,num_pix_y,aspt_ratio,'off');

%% coordinates of cytoplasma pixels
coord_in_cyto = intersect(coord_out_nuc,coord_in_cell,'rows');
% coord_in_nuc --> nucleoid pixels
% coord_out_cell --> pixels outside of cell

%% scatter plot different regions
if strcmp(fig_plot,'on')
    figure
    scatter(coord_in_nuc(:,1),coord_in_nuc(:,2),5,'filled')
    hold on
    scatter(coord_out_cell(:,1),coord_out_cell(:,2),5,'filled')
    hold on
    scatter(coord_in_cyto(:,1),coord_in_cyto(:,2),5,'filled')
    axis equal
    xlim([0 aspt_ratio*num_pix_y])
    ylim([0 num_pix_y])
    % axis off
end
>>>>>>> 5e4cf7b (Initial commit)
end