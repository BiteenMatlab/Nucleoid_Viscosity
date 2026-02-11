function Traj = BD_Simu_3D(Num, step, R, L, kb_nuc, anchor, varargin)
% 3D Brownian-dynamics simulation of a probe diffusing within effective
% space defined by 
%       • an outer rod-shaped cell inner membrane (radius R, half-length of cylinder L)
%       • an inner nucleoid shell obtained by revolving a piece-wise linear
%         functions   'x = k_i y + b_i'   between anchor points.
%
%   Traj = BD_Simu_3D(Num, step, R, L, kb_nuc, anchor, <options>)
%
% Required
%   Num        integer     – number of steps
%   step       double      – step size (same length units as R, L)
%   R, L       doubles     – cell inner membrane geometry
%   kb_nuc     N×2 double  – [k  b]  slope & intercept of each segment
%   anchor     (N+1)×2     – anchor points [y  x] (monotonically increasing y)
%
% Options (name,value, defaults shown)
%   'RandomStep'    'off' | 'on'                step size -- fixed or following rayleigh distribution      [off]
%   'Confinement'   'Collision' | 'RejectOut'   method used to introduce confinment                        [Collision]
%   'PlotTraj'      'off' | 'on'                                            [off]
%   'Gap'           positive integer – thinning factor for plotting         [200]
%   'MaxReflection' positive integer – bounce limit per collision           [25]
%   'colorMap'      color map used for trajectory visulization              [parula]
%
% Output
%   Traj  (Num+1)×7 : [x y z | step theta phi reg]
%          reg = 0 free   1 bounced at cell wall   2 bounced at nucleoid
%
% Example:
% Traj = BD_Simu_3D_v4(5e4, 0.1, 5, 6,...
%    kb_nuc_spMap, anchor_spMap, ...
%    'RandomStep','off', 'Confinement','Collision',...
%    'PlotTraj','on','colorMap',jet);
%
%
% Author: Xiaofeng Dai
% Date: 04/15/2025


%% 1 input options -----------------------------------------------------------
p = inputParser;
p.addParameter('RandomStep','off', ...
        @(s)any(strcmpi(char(s),{'on','off'})));
p.addParameter('Confinement','Collision', ...
        @(s)any(strcmpi(char(s),{'Collision','RejectOut'})));
p.addParameter('PlotTraj','off', ...
        @(s)any(strcmpi(char(s),{'on','off'})));
p.addParameter('Gap',200, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('MaxReflection',25, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('colorMap','parula', ...          
        @(c) ( (ischar(c) || isstring(c)) ...    % built-in name
              || isa(c,'function_handle') ...    % e.g. jet
              || (isnumeric(c) && size(c,2)==3) ) ); % Nx3 
p.parse(varargin{:});
opt = p.Results;
randStep = strcmpi(opt.RandomStep,'on');

%% 2 cell geometries (inner membrane and nucleoid) ---------------------------
CellX = @(y) (y>=-(R+L)&y<-L).*sqrt(-(y+L).^2+R^2) + ...
             (y>=-L & y<=L).*R + ...
             (y>L  & y<=(L+R)).*sqrt(-(y-L).^2+R^2);
NucX  = @(y) piecefunc(anchor,y);

inCell = @(p) abs(coordTrans(p)) <= CellX(p(2))+1e-10;
inNuc  = @(p) abs(coordTrans(p)) <= NucX(p(2))+1e-10;

%% 3 safe sub-step size  ---------------------------------------------------
yProbe   = linspace(-(L+R), L+R, 400);
gapAll   = CellX(yProbe) - abs(NucX(yProbe));
gap_min0 = min(gapAll(gapAll>0));           % ignore zeros
GAP_MIN  = max(step, 0.3*gap_min0);         % never smaller than one step
MAX_SUB  = 20;                              % cap sub-steps for speed

%% 4 initial coordinates ----------------------------------------------
[p0x,p0y,p0z] = randomPoint(CellX,NucX,R,L);

%% 5 allocate trajectory ----------------------------------------------
Traj        = zeros(Num+1,7);
Traj(1,1:3) = [p0x p0y p0z];

%% 6 compute displacements --------------------------------
phi   = 2*pi*rand(Num,1);
theta = pi*(rand(Num,1)-0.5);
if randStep
    alpha = raylrnd(1/sqrt(2*log(2)/log(exp(1))),Num,1);      % rayleigh distribution with median at 1  
else
    alpha = ones(Num,1);
end
len   = step*alpha;
cosT  = cos(theta);
steps = [len.*cosT.*cos(phi) , len.*cosT.*sin(phi) , len.*sin(theta)];

%% 7 introduce confinement ------------------------------------------------
for i = 1:Num
    dr   = steps(i,:);
    nSub = max(1, min(MAX_SUB, ceil(norm(dr)/GAP_MIN)));
    subD = dr / nSub;

    pos  = Traj(i,1:3);  reg = 0;
    for s = 1:nSub
        if strcmpi(opt.Confinement,'RejectOut')
            [pos,reg] = step_reject(pos, subD, inCell, inNuc, opt.MaxReflection);
        else
            [pos,reg] = step_collide(pos, subD, CellX, NucX, ...
                                     R, L, kb_nuc, anchor, opt.MaxReflection);
        end
    end
    Traj(i+1,1:3) = pos;
    Traj(i,4:7)   = [norm(dr) theta(i) phi(i) reg];
end

%% 8 plot --------------------------------------------------------------
if strcmpi(opt.PlotTraj,'on')
    plotTrajectory(Traj,R,L,NucX,anchor,opt.Gap,opt.colorMap);
end

end

%% ======================================================================
%%  Helper functions
function [x,y,z] = randomPoint(CellX,NucX,R,L)
while true
    y = -(R+L)+2*(R+L)*rand;
    x = -R     +2*R*rand;
    z = -R     +2*R*rand;
    p = [x y z];
    if abs(coordTrans(p))<=CellX(y) && abs(coordTrans(p))>NucX(y)
        return
    end
end
end
% ----------------------------------------------------------------------
function [next,reg] = step_reject(pos,dr,inCell,inNuc,maxRef)
next = pos+dr; reg = 0;
if inCell(next)&&~inNuc(next), return, end
reg = 1 + inCell(next);                     % 1 outside cell, 2 inside nuc
for t = 1:maxRef                               
    phi   = 2*pi*rand;
    theta = pi*(rand-0.5);
    dr2   = norm(dr)*[cos(theta)*cos(phi) cos(theta)*sin(phi) sin(theta)];
    cand  = pos+dr2;
    if inCell(cand)&&~inNuc(cand), next=cand; return, end
end
next = pos;                                 
end
% ----------------------------------------------------------------------
function [next,reg] = step_collide(pos,dr,CellX,NucX,R,L,kb,anchor,maxRef)
rem = dr; reg = 0;
for n = 1:maxRef
    tgt = pos + rem;
    if abs(coordTrans(tgt))<=CellX(tgt(2))+1e-10 && ...
       abs(coordTrans(tgt))> NucX(tgt(2))+1e-10
        next = tgt; return
    end
    try
        [hit,reg] = firstIntersection(pos,tgt,CellX,NucX,R,L,kb,anchor);
    catch
        next = pos; return                 % rare fallback
    end
    used = hit-pos; rem = rem-used; pos = hit;
    nrm  = surfaceNormal(hit,reg,R,L,kb,anchor);
    % rem  = rem - 2*(rem*nrm')*nrm;
    rem  = rem - 2*dot(rem,nrm)*nrm;
    if norm(rem)<1e-12, next=pos; return, end
end
next = pos;
end
% ----------------------------------------------------------------------
function [hit,reg] = firstIntersection(p0,p1,CellX,NucX,R,L,kb,anchor)
% Robust quadratic/linear solver with axial and parallel special cases.
d  = p1 - p0;
hit = nan(1,3); t_best = inf; reg = 0;
A  = d(1)^2 + d(3)^2; B = 2*(p0(1)*d(1)+p0(3)*d(3)); C = p0(1)^2+p0(3)^2-R^2;

% --- axial-only displacement ------------------------------------------
if A < 1e-14
    r2 = p0(1)^2 + p0(3)^2;
    % nucleoid crossings
    for s = 1:size(kb,1)
        k = kb(s,1); b = kb(s,2);
        rhs = sqrt(r2);
        if abs(k)<1e-14, continue, end
        y_cross = (rhs - b)/k;
        y1 = anchor(s,1); y2 = anchor(s+1,1);
        if y_cross>=y1 && y_cross<=y2
            t = (y_cross - p0(2))/d(2);
            if t>0 && t<1 && t<t_best
                hit = p0 + t*d; reg = 2; t_best = t;
            end
        end
    end
    % hemispherical caps
    if abs(p0(2))<=L && abs(p1(2))>=L
        s  = sign(d(2));                     % +1 up, -1 down
        r2 = p0(1)^2 + p0(3)^2;
        dy = sqrt(R^2 - r2);                % distance to sphere along y
        y_hit =  s*L + s*dy;                % first point on the cap
        t = (y_hit - p0(2))/d(2);
        if t>0 && t<1 && t<t_best
            hit = p0 + t*d;  reg = 1; t_best = t;
        end
    end
    if isnan(hit(1)), error('intersection fail (axial)'), end
    return
end

% --- cylinder wall -----------------------------------------------------
for t = quadraticRoots(A,B,C)'
    if 0<t && t<1
        y = p0(2)+t*d(2);
        if abs(y)<L && t<t_best
            hit = p0 + t*d; reg = 1; t_best = t;
        end
    end
end
% --- hemispherical caps ------------------------------------------------
for s = [-1 1]
    Ac = A + d(2)^2;
    Bc = B + 2*(p0(2)-s*L)*d(2);
    Cc = p0(1)^2+p0(3)^2 + (p0(2)-s*L)^2 - R^2;
    for t = quadraticRoots(Ac,Bc,Cc)'
        if 0<t && t<1 && t<t_best
            y = p0(2)+t*d(2);
            if (s==1&&y>=L)||(s==-1&&y<=-L)
                hit = p0 + t*d; reg = 1; t_best = t;
            end
        end
    end
end
% --- nucleoid segments -------------------------------------------------
for s = 1:size(kb,1)
    k = kb(s,1); b = kb(s,2); y1 = anchor(s,1); y2 = anchor(s+1,1);
    A2 = A - k^2*d(2)^2;
    B2 = B - 2*k*d(2)*(k*p0(2)+b);
    C2 = C - (k^2*p0(2)^2+2*k*b*p0(2)+b^2);
    if abs(A2)<1e-14                     % parallel → linear
        if abs(B2)<1e-14, continue, end
        t_lin = -C2/B2;
        if 0<t_lin && t_lin<1
            y = p0(2)+t_lin*d(2);
            if y>=y1&&y<=y2 && t_lin<t_best
                hit = p0 + t_lin*d; reg = 2; t_best = t_lin;
            end
        end
    else
        for t = quadraticRoots(A2,B2,C2)'
            if 0<t && t<1
                y = p0(2)+t*d(2);
                if y>=y1&&y<=y2 && t<t_best
                    hit = p0 + t*d; reg = 2; t_best = t;
                end
            end
        end
    end
end
if isnan(hit(1)), error('intersection fail (none)'), end
end
% ----------------------------------------------------------------------
function n = surfaceNormal(p,reg,R,L,kb,anchor)
if reg==1                                  % inner membrane, normal face inward
    if abs(p(2))<=L, n = -[p(1) 0 p(3)];
    else            n = -[p(1) p(2)-sign(p(2))*L p(3)];
    end
else                                        % nucleoid
    y = p(2); idx = find(y>=anchor(:,1),1,'last'); k = kb(idx,1);
    n = [p(1) -k^2*y - k*kb(idx,2) p(3)];
end
n = n./norm(n);
end
% ----------------------------------------------------------------------
function tx = coordTrans(p)
x=p(1); y=p(2); z=p(3);
if x==0 && z==0, tx = 0;
else              tx = sqrt(x.^2+z.^2).*sign(x + (x==0).*sign(z));
end
end
% ----------------------------------------------------------------------
function r = quadraticRoots(a,b,c)
disc = b.^2 - 4*a.*c; r = [];
if disc>=0
    sd = sqrt(disc);
    r  = [(-b-sd)/(2*a) ; (-b+sd)/(2*a)];
end
end
% ----------------------------------------------------------------------
function eq = piecefunc(anchor,y)
eq = zeros(size(y));
for i=1:size(anchor,1)-1
    y1=anchor(i,1); y2=anchor(i+1,1);
    k=(anchor(i+1,2)-anchor(i,2))/(y2-y1);
    b=(y2*anchor(i,2)-y1*anchor(i+1,2))/(y2-y1);
    eq = eq + ((y>=y1)&(y<y2)).*(k*y + b);
end
end
% ----------------------------------------------------------------------


