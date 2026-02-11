function [cell_outline, nuc_outline, f] = Outline_by_htMap(htMap, raw_heatmap, SmthMethod, nucMethod, aspt_ratio, num_pix_y, figplot, show_nuc, scale_factor, nuc_intp)
% Outline_by_htMap Determines cell & nucleoid boundaries from a localization heat map.
%
% [cell_outline, nuc_outline, f] = Outline_by_htMap(htMap, raw_heatmap, SmthMethod,
% nucMethod, aspt_ratio, num_pix_y, figplot, show_nuc, scale_factor, nuc_intp)
%
% This function first smooths a heat map of normalized localizations. 
% Then the gradient of the smoothed map is computed
% and used to define a baseline via histogram smoothing and clustering.
% Using this baseline, a threshold is applied to the gradient to extract a
% cell outline. Next, using a second step (and one of two nucleoid detection methods),
% a nucleoid (cytoplasm) outline is determined.
%
% Optional smoothing methods:
% 'kde' – use kernel density estimation (default)
% 'interpolation' – use cubic interpolation followed by Gaussian smoothing.
%
% Nucleoid detection methods:
% 'gradient' – uses local gradient information (default)
% 'threshold' – uses a threshold based on the midline intensity.
%
% INPUTS (with defaults if not provided):
% htMap - Original heat map (matrix) before additional smoothing.
% raw_heatmap- [N x 2] matrix of normalized localization coordinates.
% SmthMethod - (default 'kde') Method for smoothing the raw data.
% nucMethod - (default 'gradient') Method for nucleoid detection.
% aspt_ratio - (default 2) Cell aspect ratio.
% num_pix_y - (default 20) Number of grid bins in y–direction.
% figplot - (default 'on') 'on' plots the resulting maps and outlines.
% show_nuc - (default 'on') 'on' overlays the nucleoid outline.
% scale_factor - (default 10) Scaling factor for interpolation resolution.
% nuc_intp - (default 'off') 'on' forces interpolation of nucleoid outline.
%
% OUTPUTS:
% cell_outline - [M x 2] coordinates (normalized) outlining the cell.
% nuc_outline - [K x 2] coordinates outlining the nucleoid/cytoplasm region.
% f - (Only used in the 'gradient' nucleoid method) Additional diagnostic
% data computed for nucleoid detection.
%
% Author: Xiaofeng Dai
% Date: 04/12/2025

%% 1. Set Default Parameters (if not provided)
if nargin < 3 || isempty(SmthMethod), SmthMethod = 'kde'; end
if nargin < 4 || isempty(nucMethod),  nucMethod  = 'gradient'; end
if nargin < 5 || isempty(aspt_ratio), aspt_ratio = 2; end
if nargin < 6 || isempty(num_pix_y),  num_pix_y  = 20; end
if nargin < 7 || isempty(figplot),    figplot    = 'on'; end
if nargin < 8 || isempty(show_nuc),     show_nuc   = 'on'; end
if nargin < 9 || isempty(scale_factor), scale_factor = 10; end
if nargin < 10 || isempty(nuc_intp),    nuc_intp   = 'off'; end

%% 2. Smooth the Heat Map
% Create a smoothed map (map_smth) using the selected smoothing method.
switch lower(SmthMethod)
    case 'kde'
        % Replicate raw data in four quadrants to force symmetry.
        data1 = raw_heatmap(:, 1:2);
        data2 = -data1;
        data3 = data1; data3(:,1) = -data1(:,1);
        data4 = -data3;
        data_loc = cat(1, data1, data2, data3, data4);
        % Define grid for density estimation.
        gridx = linspace(-0.5, 0.5, num_pix_y * aspt_ratio * scale_factor);
        gridy = linspace(-0.5, 0.5, num_pix_y * scale_factor);
        [X, Y] = meshgrid(gridx, gridy);
        pts = [X(:), Y(:)];
        % Estimate density.
        [map_kde,~] = ksdensity(data_loc, pts);
        map_kde = reshape(map_kde, [num_pix_y * scale_factor, num_pix_y * aspt_ratio * scale_factor]);
        map_smth = map_kde;
    case 'interpolation'
        % Interpolate the original heat map and smooth it.
        map_intp = Map_Interpolate(htMap, scale_factor, 'cubic');
        map_intp_smth = smoothdata2(map_intp,"gaussian",SmoothingFactor=0.03);
        map_smth = map_intp_smth;
    otherwise
        error('Unknown smoothing method: %s', SmthMethod);
end

%% 3. Compute Gradient and Determine Baseline via Histogram Analysis
[g1, g2] = gradient(map_smth);
G = sqrt(g1.^2 + g2.^2);
% Flatten gradient image into a vector.
J = G(:);
J(isnan(J)) = [];
% Compute histogram of gradient magnitudes.
[occ_his, occ_edges] = histcounts(J, 75, 'Normalization', 'probability');
occ_bin = Get_bincenter(occ_edges);
% Remove first bin (often an outlier).
occ_his(1) = []; occ_bin(1) = [];
% Smooth the occurrence histogram.
occ_his_smth = smooth(occ_his, 8, 'moving');
occ_grd = smooth(gradient(occ_his_smth));
% Find indices where the gradient of the occurrence histogram falls within a narrow range.
occ_flat = find(occ_grd > 0.15 * min(occ_grd) & occ_grd < 0.15 * max(occ_grd));
buff = zeros(length(occ_flat), 2);
buff(:,1) = occ_flat;
occ_flat = buff; clear buff
% Cluster these indices using DBSCAN.
occ_flat_ind = dbscan(occ_flat, 1, 2);
occ_flat(:,2) = occ_flat_ind;
occ_flat(occ_flat(:,2)==-1, :) = [];  % remove noise.
sel_group = mode(occ_flat(:,2));       % select group with most dots contained
base_level = mean(occ_his_smth(occ_flat(occ_flat(:,2)==sel_group, 1)));
occ_his_smth = occ_his_smth - base_level;
occ_his_smth(occ_his_smth < 0) = 0;
% Fit the smoothed histogram with a dual–Gaussian model.
para = G2P(occ_bin, occ_his_smth);
% Set a threshold value for the gradient image.
thres_val = para(2) + 2*para(3);          % para(2) + 1.75*para(3); para(2) + 1*para(3)
clear occ_edges occ_his occ_bin occ_flat occ_grd base_level

%% 4. Extract Cell Outline from Gradient Image
H = G;
H(H < thres_val) = 0;
H(isnan(H)) = 0;
% Get coordinates of nonzero gradient points.
[coordy, coordx] = find(H ~= 0);
coord = [coordx, coordy];
% Cluster the coordinates to select the main group.
coord_clus = dbscan(coord, 1, 5);
coord(coord_clus < 1, :) = [];
coord_clus(coord_clus<1) = [];

bndy = coord(coord_clus == 1, :);
% Compute a boundary (with shrink factor 0.2) and convert to normalized coordinates.
k = boundary(bndy(:,1), bndy(:,2), 0.2);
cell_outline = bndy(k, :) ./ scale_factor + 0.5;
clear coord coord_clus k H
%% 5. Determine Cytoplasm / Nucleoid Outline
k = boundary(bndy(:,1),bndy(:,2));
bndy1 = bndy(k,:);
H_reg = InRegion(G,bndy1);

thres_val2 = para(5) - para(6);
H_reg(H_reg > thres_val2) = 0;
H_reg(isnan(H_reg)) = 0;
[test1, test2] = find(H_reg ~= 0);
test_coords = [test1, test2];
k = dbscan(test_coords, 2, 8);
clear test1 test2
gp = mode(k);
buff = test_coords(k == gp, :);
% Compute a boundary for the nucleoid region.
bndy2 = boundary(buff(:,1), buff(:,2));
bndy2 = buff(bndy2, :);
% Swap columns if needed (to get proper [x,y] order)
tmp = bndy2(:,1);
bndy2(:,1) = bndy2(:,2);
bndy2(:,2) = tmp;
clear tmp k

% Based on the chosen nucleoid detection method…
switch lower(nucMethod)
    case 'threshold'
        % Use midline intensity to define a threshold.
        mid_int = map_smth(:,round(4*width(map_smth)/8));
        % Find the index corresponding to the maximum intensity.
        [~, idx_max] = max(mid_int);
        c = length(mid_int) - idx_max + 1;
        lb = min(idx_max, c);
        rb = max(idx_max, c);
        mid_int = mid_int(lb:rb);
        mid_g = gradient(mid_int);
        % Define the nucleoid threshold as the average intensity at the minimum and maximum gradient.
        nuc_thresh = (min(mid_int(mid_g==min(mid_g))) + max(mid_int(mid_g==max(mid_g)))) / 2;
        Nuc_region = InRegion(map_smth, bndy2);
        Nuc_region(map_smth >= nuc_thresh) = 0;
        Nbndy = GetRegion(Nuc_region);
        k = boundary(Nbndy(:,1), Nbndy(:,2), 0.25);
        nuc_outline = Nbndy(k, :) ./ scale_factor + 0.5;
        clear bndy1 H_reg mid_g mid_int Nuc_region k nuc_thresh lb rb c
    case 'gradient'
        % Use the gradient along the midline to estimate a boundary.
        cnt_int = map_smth(round(height(map_smth)/2), :);
        [xs, xe, ~] = GRange(cnt_int, 'FWHM');
        Ind = xs:xe;
        d = zeros(length(Ind), 2);
        e = zeros(2 * length(Ind), 2);
        for i = 1:length(Ind)
            ind = Ind(i);
            col_data = map_smth(:, ind);
            [a, b, ~] = GRange(col_data, 'MaxGradient');
            d(i,:) = [a, b];
            e(2*i-1,:) = [ind, a];
            e(2*i,:)   = [ind, b];
        end
        ys = max(d(:,1));
        ye = min(d(:,2));
        clear d
        ctplsm = InRegion(map_smth, bndy2);
        [g1, g2] = gradient(ctplsm);
        G_vec = sqrt(g1.^2 + g2.^2);
        clear g1 g2 ctplsm
        J = G_vec(:);
        J(isnan(J)) = [];
        J(J==0) = [];
        outliers_removed = rmoutliers(J);
        G_vec(G_vec > max(outliers_removed)) = 0;
        clear J
        G_smth = smoothdata2(G_vec, "gaussian", 'SmoothingFactor', 0.08);
        [occ_his, occ_bin] = histcounts(outliers_removed, 75, 'Normalization', 'probability'); 
        occ_bin = Get_bincenter(occ_bin);
        occ_his(1) = [];
        occ_bin(1) = [];
        occ_his_smth = smooth(occ_his, 4, 'moving');
        occ_grd = smooth(gradient(occ_his_smth));
        occ_flat = find(occ_grd > 0.15 * min(occ_grd) & occ_grd < 0.15 * max(occ_grd));
        buff = zeros(length(occ_flat),2);
        buff(:,1) = occ_flat;
        occ_flat = buff; clear buff
        occ_flat_ind = dbscan(occ_flat, 1, 2);
        occ_flat(:,2) = occ_flat_ind;
        occ_flat(occ_flat(:,2)==-1, :) = [];
        sel_group = mode(occ_flat(:,2));
        thres = occ_bin(occ_flat(find(occ_flat(:,2)==sel_group,1,'first'),1) - min(1, round(length(occ_bin)*0.05)));
        gg = G_smth;
        gg(gg > thres) = 0;
        f = zeros(2 * length(ys:ye), 2);
        cnt = 1;
        for i = ys:ye
            row_data = smoothdata(gg(i, :), 'gaussian', round(length(gg(i, :))*0.03));
            a_val = find(row_data == max(row_data), 1, 'first');
            b_val = num_pix_y * aspt_ratio * scale_factor - a_val;
            f(2*cnt-1,:) = [min(a_val, b_val), i];
            f(2*cnt, :)  = [max(a_val, b_val), i];
            cnt = cnt + 1;
        end
        buff = f(2:2:size(f,1), :);
        idx = dbscan(buff, 2, 2);
        bb = [num_pix_y * aspt_ratio * scale_factor - min(buff(idx == mode(idx),1)), ...
              min(buff(idx == mode(idx), 1))];
        % Filter out any points outside the desired column range.
        e(e(:,1) < bb(1) | e(:,1) > bb(2), :) = [];
        ef = cat(1, e, f);
        k = boundary(ef(:,1), ef(:,2), 0.25);
        Nbndy = ef(k, :);
        nuc_outline = Nbndy ./ scale_factor + 0.5;
    otherwise
        error('Wrong Nucleoid Detection Method: %s', nucMethod);
end

%% 6. Optionally Interpolate the Nucleoid Outline
if strcmpi(nuc_intp, 'on')
    nuc_outline = Intp_Nuc(nuc_outline, num_pix_y, 10);
end

%% 7. Display Figure with Overlaid Outlines (if requested)
if strcmpi(figplot, 'on')
    f = figure;
    imagesc(htMap);
    if strcmpi(show_nuc, 'on')
        hold on;
        if strcmpi(nuc_intp, 'off')
            plot(nuc_outline(:,1), nuc_outline(:,2), 'red', 'LineWidth', 1.2, 'LineStyle', '--');
            % Connect first and last point to complete the outline.
            plot([nuc_outline(end,1); nuc_outline(1,1)], [nuc_outline(end,2); nuc_outline(1,2)], 'red');
        else
            scatter(nuc_outline(1:5:end,1), nuc_outline(1:5:end,2), 1.6, 'red', 'filled');
        end
    end
    hold on;
    % Plot cell outline in white.
    plot(cell_outline(:,1), cell_outline(:,2), 'white', 'LineWidth', 1.5);
    hold on;
    plot([cell_outline(end,1); cell_outline(1,1)], [cell_outline(end,2); cell_outline(1,2)], 'white', 'LineWidth', 1.5);
    pbaspect([aspt_ratio, 1, 1]);
    num_pix_x = num_pix_y * aspt_ratio;
    xlim([0.5, num_pix_x + 0.5]);
    xticks(linspace(0.5, num_pix_x + 0.5, 6));
    ylim([0.5, num_pix_y + 0.5]);
    yticks(linspace(0.5, num_pix_y + 0.5, 6));
    xticklabels({'0','0.2','0.4','0.6','0.8','1'});
    yticklabels({'1','0.8','0.6','0.4','0.2','0'});
    colormap(slanCM('viridis'))
    colorbar;
else
    f = [];  % No figure handle if not plotted.
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunction: Get_bincenter
function Bincenter = Get_bincenter(edges)
% Computes bin center positions from histogram bin edges.
Bincenter = (edges(1:end-1) + edges(2:end)) / 2;
Bincenter = Bincenter.';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunction: InRegion
function in_region = InRegion(G, bndy)
% Returns the values of matrix G that fall within the polygon defined by bndy.
% bndy is an [N x 2] matrix of (x,y) coordinates.
mask = poly2mask(bndy(:,1), bndy(:,2), size(G,1), size(G,2));
in_region = G .* mask;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunction: GetRegion
function region = GetRegion(H)
% Extracts the main region (cluster) from nonzero entries in H.
[r, c] = find(H ~= 0);
coord = [c, r];
clusters = dbscan(coord, 2, 8);
% Remove noise and select the most frequent cluster.
coord(clusters < 1, :) = [];
clusters(clusters < 1) = [];
sel_group = mode(clusters);
region = coord(clusters == sel_group, :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunction: GRange

function [thres1, thres2, intensity] = GRange(dis,method)
if nargin==1
    method = 'MaxGradient';
end
b = find(dis==max(dis),1,"first");
c = length(dis)-b+1;
lb = min(b,c);
rb = max(b,c);
clear b c
dis = dis(lb:rb);
clear rb
switch method
    case 'FWHM'
        buff = (max(dis) + min(dis))/2;
        [~, ind] = min(abs(dis-buff));
        clear buff
        b = ind;
        clear ind
        c = length(dis)-b+1;
        thres1 = min(b,c)+lb;
        thres2 = max(b,c)+lb;
        intensity = (dis(b)+dis(c))/2;
        clear b c
    case 'MaxGradient'
        dis_g = gradient(dis);
        thres1 = find(dis_g==min(dis_g))+lb;
        thres2 = find(dis_g==max(dis_g))+lb;
        intensity = (dis(dis_g==min(dis_g))+dis(dis_g==max(dis_g)))/2;
    otherwise
        error('Wrong Method')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunction: G2P
function [p_fit, err] = G2P(xdata, ydata)
% G2P fits a two-peak Gaussian model to the data.
% Model: p(1)*exp(-((x-p(2))/p(3)).^2) + p(4)*exp(-((x-p(5))/p(6)).^2) + p(7)
fun = @(p,xdata)  p(1)*exp(-((xdata - p(2))/p(3)).^2) + p(4)*exp(-((xdata - p(5))/p(6)).^2) + p(7);
p_initial = [0.01, (mean(xdata)+min(xdata))/2, mean(xdata)/20, ...
0.01, (mean(xdata)+max(xdata))/2, max(xdata)/10, min(ydata)/20];
lb = [0, min(xdata), 0, 0, min(xdata), 0, 0];
ub = [1, max(xdata), max(xdata), 1, max(xdata), max(xdata), 0.1*max(ydata)];
options = optimset('Display','off');
[p_fit, err] = lsqcurvefit(fun, p_initial, xdata, ydata, lb, ub, options);
end