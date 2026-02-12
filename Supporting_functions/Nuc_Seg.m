function [coor_in_core, coor_in_peri, peri1, peri2] =  Nuc_Seg(spMap, contour_range, cell_outline, num_pix_y, aspt_ratio, figPlot)

if nargin < 4 || isempty(num_pix_y),    num_pix_y  = 20;             end
if nargin < 5 || isempty(aspt_ratio),   aspt_ratio = 2;              end
if nargin < 6 || isempty(figPlot),      figPlot    = ["off" "off"];  end

new_nuc_contour = Outline_by_spMap(spMap,contour_range,...
    figPlot(1),cell_outline,aspt_ratio,num_pix_y);
peri1 = new_nuc_contour{1};
peri2 = new_nuc_contour{2};

[wcoord,lcoord] =...
    meshgrid((1:1:height(spMap)),(1:1:width(spMap)));
map_coord = cat(2,lcoord(:),wcoord(:));

ind_peri = inpolygon(map_coord(:,1),map_coord(:,2),peri1(:,1),peri1(:,2));
ind_core = inpolygon(map_coord(:,1),map_coord(:,2),peri2(:,1),peri2(:,2));
a = find(ind_core==0);
b = find(ind_peri==1);
c = intersect(a,b);
clear a b
coor_in_core = map_coord(ind_core==1,:);
coor_in_peri = map_coord(c,:);
clear c

if strcmpi(figPlot(2),'on')
    figure
    scatter(coor_in_peri(:,1),coor_in_peri(:,2),3,'blue','filled')
    hold on
    scatter(coor_in_core(:,1),coor_in_core(:,2),6,'red','filled')
    hold on
    plot(polyshape(peri1(:,1),peri1(:,2)),'FaceColor','#D95319','EdgeAlpha',0.05,'EdgeColor','#D95319')
    hold on
    plot(polyshape(peri2(:,1),peri2(:,2)),'FaceColor','#D95319','EdgeAlpha',0.05,'EdgeColor','#D95319')
    hold on
    plot(cell_outline(:,1), cell_outline(:,2), 'black', 'LineWidth', 0.75);
    % Connect first and last point.
    plot([cell_outline(end,1); cell_outline(1,1)], [cell_outline(end,2); cell_outline(1,2)], 'black', 'LineWidth', 0.75);
    axis equal
    xlim([0 aspt_ratio*num_pix_y])
    ylim([0 num_pix_y])
end

end