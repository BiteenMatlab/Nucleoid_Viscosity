function plotTrajectory(Traj,R,L,NucX,anchor,gap,cMap)
keep = 1:gap:size(Traj,1); 
% cmap = parula(numel(keep));
if isa(cMap,'function_handle')         % e.g. @turbo
    cmap = cMap(numel(keep));
elseif ischar(cMap) || isstring(cMap)  % e.g. 'jet'
    cmap = feval(char(cMap),numel(keep));
else                                   % numeric Nx3 matrix passed in
    cmap = cMap;
    % if its length differs from keep, interpolate:
    if size(cmap,1)~=numel(keep)
        cmap = interp1( linspace(0,1,size(cmap,1)), ...
                        cmap, linspace(0,1,numel(keep)) );
    end
end
figure('Color','w'); hold on; view(3); grid on; axis equal
xlabel('X', 'FontWeight', 'bold');
ylabel('Y', 'FontWeight', 'bold');
zlabel('Z', 'FontWeight', 'bold');
% xlabel X; ylabel Y; zlabel Z; 
title('3D Brownian Dynamics Simulation');
% ----- cell cylinder ---------------------------------------------------
n = 100; [th,yc] = meshgrid(linspace(0,2*pi,n), linspace(-L,L,n));
Xc = R*cos(th); Zc = R*sin(th); Yc = yc;
surf(Xc,Yc,Zc,'EdgeColor','none','FaceAlpha',0.1,'FaceColor',[.5 .5 .5]);
% ----- hemispherical caps ---------------------------------------------
[thC,phi] = meshgrid(linspace(0,2*pi,n), linspace(0,pi/2,n/2));
Xcap = R*sin(phi).*cos(thC); Zcap = R*sin(phi).*sin(thC);
YcapTop = L + R*cos(phi); YcapBot = -L - R*cos(phi);
surf(Xcap,YcapTop,Zcap,'EdgeColor','none','FaceAlpha',0.1,'FaceColor',[.5 .5 .5]);
surf(Xcap,YcapBot,Zcap,'EdgeColor','none','FaceAlpha',0.1,'FaceColor',[.5 .5 .5]);
% ----- nucleoid shell --------------------------------------------------
Sn = nucSurface(NucX,anchor,n,n);
surf(Sn(:,:,1),Sn(:,:,2),Sn(:,:,3), ...
     'EdgeColor','none','FaceAlpha',0.35,'FaceColor',[0 .45 .74]);
% ----- trajectory ------------------------------------------------------
for k = 1:numel(keep)-1
    plot3(Traj(keep(k:k+1),1), Traj(keep(k:k+1),2), Traj(keep(k:k+1),3), ...
          'LineWidth',1,'Color',cmap(k,:));
end
colormap(cmap); colorbar;

% ----- axises set up ---------------------------------------------------
scalingCellWidth = R*1.05;
scalingCellLength = (R+L)*1.05;
xlim([-scalingCellWidth, scalingCellWidth]);
ylim([-scalingCellLength, scalingCellLength]);
zlim([-scalingCellWidth, scalingCellWidth]);
xticks([-scalingCellWidth, 0, scalingCellWidth]);
xticklabels({'-0.5','0','0.5'});
yticks([-scalingCellLength,-0.5*scalingCellLength, 0,0.5*scalingCellLength, scalingCellLength]);
yticklabels({'-0.5','-0.25','0','0.25','0.5'});
zticks([-scalingCellWidth, 0, scalingCellWidth]);
zticklabels({'-0.5','0','0.5'});

nSeg  = numel(keep);           % # line segments (colour indices go 1‥nSeg)
clim([1 nSeg]);               % map colour index → axis on the bar
cb    = colorbar;              % handle comes from last call to colorbar()
stepVals   = 0 : 75*gap : (size(Traj,1)-1);
tickIdx    = stepVals/gap + 1;            % convert to colour-index positions
valid      = tickIdx <= nSeg;             % don’t exceed the bar’s range
cb.Ticks   = tickIdx(valid);              % tick positions (in colour index units)
cb.TickLabels = arrayfun(@num2str, stepVals(valid)/gap, 'UniformOutput',false);

end
% ----------------------------------------------------------------------
function Sn = nucSurface(NucX,anchor,ny,nt)
y = linspace(anchor(1,1), anchor(end,1), ny);
x_half = arrayfun(NucX,y);
th = linspace(0,2*pi,nt)';
[Y,TH] = meshgrid(y,th); [Xh,~] = meshgrid(x_half,th);
Sn(:,:,1) = Xh.*cos(TH);
Sn(:,:,3) = Xh.*sin(TH);
Sn(:,:,2) = Y;
end