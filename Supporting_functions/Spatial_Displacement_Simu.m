function [Simu_spMap,Simu_spatialDisp] = Spatial_Displacement_Simu(Store_result,gap,num_pix_y,aspt_ratio,scaling_width,scaling_length,scale_bar,thres_ratio,symmetrize,fig_spMap)
if nargin==7
    thres_ratio = 0.05;
    symmetrize = 'on';
    fig_spMap = slanCM('torch');
end


Num_run = length(Store_result);
Traj_cat = cell(Num_run,1);
parfor j = 1:Num_run
    traj = Store_result{j};
    traj_long = [traj(1:gap:height(traj),1),traj(1:gap:height(traj),2),traj(1:gap:height(traj),3)];
    traj_project = cat(2,(1:1:height(traj_long)).',traj_long(:,1:2));   %projection on x-y plane, with frame index col1
    Traj_cat{j} = traj_project;
end
clear i traj traj_long traj_project

num_pix_x = num_pix_y*aspt_ratio;
NormSSS = SpatialStepSize(cat(1,Traj_cat{:}));     % calculate Spatial Step Size 
% normalization
NormSSS(:,1) = NormSSS(:,1)/scaling_width/2;
NormSSS(:,2) = NormSSS(:,2)/scaling_length/2;
NormSSS(:,3) = NormSSS(:,3)*scale_bar;
buff_w = NormSSS(:,1);
buff_l = NormSSS(:,2);
NormSSS(:,1) = buff_l;
NormSSS(:,2) = buff_w;
clear buff_l buff_w
Simu_spatialDisp = NormSSS;

if strcmp(symmetrize,'off') % Unsymetrized
    [DisZ,CntZ]=PlotSSS(NormSSS,num_pix_y,aspt_ratio);
    aa = reshape(CntZ,[],1);
    bb = maxk(aa,round(0.1*length(aa)));
    threshold=thres_ratio*mean(bb);
    % threshold=thres_ratio*max(max(CntZ));
    index=CntZ<threshold;
    DisZ(index)=NaN;
    Simu_spMap = DisZ;
    if strcmp(fig_spMap,'off')==false
        figure
        imagesc(Simu_spMap);
        axis equal
        colormap(gca,fig_spMap);
        colorbar;
        xlim([0.5 num_pix_x+0.5])
        xticks([0.5 (num_pix_x/5)+0.5 2*(num_pix_x/5)+0.5 3*(num_pix_x/5)+0.5 4*(num_pix_x/5)+0.5 num_pix_x+0.5])
        ylim([0.5 num_pix_y+0.5])
        yticks([0.5 (num_pix_y/5)+0.5 2*(num_pix_y/5)+0.5 3*(num_pix_y/5)+0.5 4*(num_pix_y/5)+0.5 num_pix_y+0.5])
        xticklabels({'0','0.2','0.4','0.6','0.8','1'});
        yticklabels({'1','0.8','0.6','0.4','0.2','0'})
        pbaspect([aspt_ratio 1 1])
        title('Spatial Displacement Map from Simulation')
    end
elseif strcmp(symmetrize,'on') % Symetrized
    [DisZ,CntZ]=PlotSSS(NormSSS,num_pix_y,aspt_ratio*2);

    Sum_DisZ=DisZ.*CntZ;
    Map_Shape=size(DisZ);
    Region1=Sum_DisZ(1:(Map_Shape(1)/2),1:(Map_Shape(2)/2)); %upper-left
    Count1=CntZ(1:(Map_Shape(1)/2),1:(Map_Shape(2)/2));
    Region2=Sum_DisZ((Map_Shape(1)/2+1):Map_Shape(1),1:(Map_Shape(2)/2)); %bottom-left
    Count2=CntZ((Map_Shape(1)/2+1):Map_Shape(1),1:(Map_Shape(2)/2));
    Region2=flipud(Region2);
    Count2=flipud(Count2);
    Region3=Sum_DisZ(1:(Map_Shape(1)/2),(Map_Shape(2)/2+1):Map_Shape(2)); %upper-right
    Count3=CntZ(1:(Map_Shape(1)/2),(Map_Shape(2)/2+1):Map_Shape(2));
    Region3=fliplr(Region3);
    Count3=fliplr(Count3);
    Region4=Sum_DisZ((Map_Shape(1)/2+1):Map_Shape(1),(Map_Shape(2)/2+1):Map_Shape(2)); %bottom-right
    Count4=CntZ((Map_Shape(1)/2+1):Map_Shape(1),(Map_Shape(2)/2+1):Map_Shape(2));
    Region4=rot90(Region4,2);
    Count4=rot90(Count4,2);
    Regions=Region1+Region2+Region3+Region4;
    Counts=Count1+Count2+Count3+Count4;
    Map1=Regions./Counts;
    Map2=flipud(Map1);
    Map3=fliplr(Map1);
    Map4=rot90(Map1,2);
    Map=[Map1,Map3;Map2,Map4];
    Count_Map=[Counts,fliplr(Counts);flipud(Counts),rot90(Counts,2)];
    
    buff_map = zeros(height(Map),width(Map)/2);
    buff_count = zeros(height(Map),width(Map)/2);
    for cc = 1:width(Map)/2
        map_1 = Map(:,2*cc-1);
        map_2 = Map(:,2*cc);
        cnt_1 = Count_Map(:,2*cc-1);
        cnt_2 = Count_Map(:,2*cc);
        buff_map(:,cc) = (map_1.*cnt_1 + map_2.*cnt_2)./(cnt_1+cnt_2);
        buff_count(:,cc) = (cnt_1+cnt_2);
    end
    Count_Map = buff_count;
    Simu_spMap = buff_map;
    aa = reshape(Count_Map,[],1);
    bb = maxk(aa,round(0.1*length(aa)));
    threshold=thres_ratio*mean(bb);
    % threshold=thres_ratio*max(max(Count_Map));
    index=Count_Map<threshold;
    Simu_spMap(index)=NaN;

    if strcmp(fig_spMap,'off')==false
        figure
        imagesc(Simu_spMap);
        axis equal
        colormap(gca,fig_spMap);
        % colormap(gca,'hot');
        colorbar;
        
        xlim([0.5 num_pix_x+0.5])
        xticks([0.5 (num_pix_x/5)+0.5 2*(num_pix_x/5)+0.5 3*(num_pix_x/5)+0.5 4*(num_pix_x/5)+0.5 num_pix_x+0.5])
        ylim([0.5 num_pix_y+0.5])
        yticks([0.5 (num_pix_y/5)+0.5 2*(num_pix_y/5)+0.5 3*(num_pix_y/5)+0.5 4*(num_pix_y/5)+0.5 num_pix_y+0.5])
        xticklabels({'0','0.2','0.4','0.6','0.8','1'});
        yticklabels({'1','0.8','0.6','0.4','0.2','0'})
        pbaspect([aspt_ratio 1 1])
        title('Spatial Displacement Map from Simulation')
    end

end

end

function [SSS]=SpatialStepSize(trackinfo)
    numrow=size(trackinfo,1);
    SSS=zeros(numrow-1,3);
    for i=1:numrow-1
        if trackinfo(i,1)-trackinfo(i+1,1)==-1
            coord_row=0.5*(trackinfo(i,2)+trackinfo(i+1,2));
            coord_col=0.5*(trackinfo(i,3)+trackinfo(i+1,3));
            stepsize=((trackinfo(i,2)-trackinfo(i+1,2))^2+(trackinfo(i,3)-trackinfo(i+1,3))^2)^0.5;
            SSS(i,1)=coord_row;
            SSS(i,2)=coord_col;
            SSS(i,3)=stepsize;
        end
    end
    test = sum(SSS,2);
    ind = test==0;
    SSS(ind,:) = [];
end

% Plot Spatial StepSize distribution 
function [DisZ,CntZ]=PlotSSS(NormSSS,n,aspt_ratio)
    xscale=1/(aspt_ratio*n);
    yscale=1/n;
    DisZ=zeros(n,aspt_ratio*n);
    CntZ=zeros(n,aspt_ratio*n);
    num=size(NormSSS,1);
    for i=1:num
        xindex=floor((NormSSS(i,1)+0.5)/xscale)+1;
        yindex=floor((NormSSS(i,2)+0.5)/yscale)+1;
        if xindex>0 && yindex>0 && xindex<= 1/xscale && yindex<= 1/yscale
            DisZ(yindex,xindex)=DisZ(yindex,xindex)+NormSSS(i,3);
            CntZ(yindex,xindex)=CntZ(yindex,xindex)+1;
        end
    end
    DisZ=DisZ./CntZ;
end