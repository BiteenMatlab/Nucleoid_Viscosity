<<<<<<< HEAD
function Aug_Result = SimuTraj_RotAug(Store_result,num_agl_seg)
if nargin==1
    num_agl_seg = 6;
end
rot_agl_rng = (1:1:num_agl_seg-1)*pi/num_agl_seg;
for i = 1:size(Store_result,2)
    Traj_sel = Store_result{1,i};
    parfor k = 1:num_agl_seg-1
        Traj_rot = zeros(size(Traj_sel));
        rot_agl = rot_agl_rng(k);
        for j = 1:size(Traj_sel,1)
            x = Traj_sel(j,1);
            y = Traj_sel(j,2);
            z = Traj_sel(j,3);
            [a1, a2, a3] = CoordRotate(x,y,z,rot_agl);
            Traj_rot(j,:) = [a1, a2, a3];
        end
        Store_result{k+1,i} = Traj_rot;
    end
end
Aug_Result = reshape(Store_result,[],1);
=======
function Aug_Result = SimuTraj_RotAug(Store_result,num_agl_seg)
if nargin==1
    num_agl_seg = 6;
end
rot_agl_rng = (1:1:num_agl_seg-1)*pi/num_agl_seg;
for i = 1:size(Store_result,2)
    Traj_sel = Store_result{1,i};
    parfor k = 1:num_agl_seg-1
        Traj_rot = zeros(size(Traj_sel));
        rot_agl = rot_agl_rng(k);
        for j = 1:size(Traj_sel,1)
            x = Traj_sel(j,1);
            y = Traj_sel(j,2);
            z = Traj_sel(j,3);
            [a1, a2, a3] = CoordRotate(x,y,z,rot_agl);
            Traj_rot(j,:) = [a1, a2, a3];
        end
        Store_result{k+1,i} = Traj_rot;
    end
end
Aug_Result = reshape(Store_result,[],1);
>>>>>>> 5e4cf7b (Initial commit)
end