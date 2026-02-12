function [coord_in_region, coor_out_region] = MapCoord_region(Map,bndy_coord,num_pix_y,aspt_ratio,fig_plot)
%% coordinates of each pixel centered at (0,0)
[wcoord,lcoord] = meshgrid((1:1:height(Map)),(1:1:width(Map)));
wcoord = wcoord-num_pix_y/2-0.5;
lcoord = lcoord-num_pix_y*aspt_ratio/2-0.5;
map_coord = cat(2,lcoord(:),wcoord(:));
map_coord(map_coord(:,2)<0,:)=[];
clear wcoord lcoord

%% coordinates in polygon defined by boundary dots
test = inpolygon(map_coord(:,1),map_coord(:,2),bndy_coord(:,1),bndy_coord(:,2));
coord_in_region = map_coord(test==1,:);
buff = coord_in_region;
% flip coordinates (symetry)
buff(:,2) = -buff(:,2);
coord_in_region = cat(1,coord_in_region,buff);
clear buff

%% coordinates outside of polygon defined by boundary dots
coor_out_region = map_coord(test==0,:);
buff = coor_out_region;
buff(:,2) = -buff(:,2);
coor_out_region = cat(1,coor_out_region,buff);
clear buff test

%% scatter plot in/out regions
if strcmp(fig_plot,'on')
    figure
    scatter(coor_out_region(:,1),coor_out_region(:,2),3,'blue','filled')
    hold on
    scatter(coord_in_region(:,1),coord_in_region(:,2),6,'red','filled')
    hold on
    warning('off','all')
    b = polyshape(bndy_coord(:,1),bndy_coord(:,2));
    warning('off','all')
    b2 = polyshape(bndy_coord(:,1),-bndy_coord(:,2));
    plot(b,'FaceColor','#D95319','EdgeAlpha',0.05,'EdgeColor','#D95319')
    hold on
    plot(b2,'FaceColor','#D95319','EdgeAlpha',0.05,'EdgeColor','#D95319')
    axis equal
    clear b b2

end

%% Coordinates start from (0,0)
coord_in_region(:,1) = coord_in_region(:,1)+num_pix_y*aspt_ratio/2+0.5;
coord_in_region(:,2) = coord_in_region(:,2)+num_pix_y/2+0.5;
coor_out_region(:,1) = coor_out_region(:,1)+num_pix_y*aspt_ratio/2+0.5;
coor_out_region(:,2) = coor_out_region(:,2)+num_pix_y/2+0.5;

end