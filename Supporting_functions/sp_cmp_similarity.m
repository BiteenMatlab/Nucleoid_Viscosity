function [sp_L2_unweighted, sp_L2, sp_L2_adj, sp_SSIM, sp_PCC, sp_JSD] = sp_cmp_similarity(spMap, Simu_spMap, weight_matrix, coord_in_nuc)
% unweighted L2 norm distance
sp_L2_unweighted = sum(sqrt((spMap - Simu_spMap).^2),'all','omitnan')/sum(~isnan(sqrt((spMap - Simu_spMap).^2)),'all');

% weighted L2 norme distance 
sp_L2 = sum(sqrt((spMap - Simu_spMap).^2).*weight_matrix,'all','omitnan');

% weighted L2 norme distance without cell pole
al = min(coord_in_nuc(:,1));
ar = max(coord_in_nuc(:,1));
spMap_adj = spMap(:,al:ar);
Simu_spMap_adj = Simu_spMap(:,al:ar);
weight_matrix_adj = weight_matrix(:,al:ar);
weight_matrix_adj = weight_matrix_adj/sum(weight_matrix_adj,"all");
sp_L2_adj = sum(sqrt((spMap_adj - Simu_spMap_adj).^2).*weight_matrix_adj,'all','omitnan');

clear al ar spMap_adj Simu_spMap_adj weight_matrix_adj

% weighted PCC
sp_PCC = weightedPCC(spMap,Simu_spMap,weight_matrix);
% weighted SSIM
sp_SSIM = weightedSSIM(spMap,Simu_spMap,weight_matrix);
% weighted JSD
sp_JSD = log(weightedJSD(spMap,Simu_spMap,weight_matrix));

end

