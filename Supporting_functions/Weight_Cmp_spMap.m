function weight_matrix = Weight_Cmp_spMap(spMap_cnt,Region,weight1,weight2,method)
if nargin == 2
    weight1 = [0.5, 1, 0];
    weight2 = [1.2, 1, 1];
    method = 'linear';
end

if nargin == 3
    weight2 = [1.2, 1, 1];
    method = 'linear';
end

if nargin == 4
    method = 'combine';
end

coord_in_nuc = Region{1};
coord_in_cyto = Region{2};
coord_out_cell = Region{3};

% normalize number of steps from each coordinate
weight_by_cnt = spMap_cnt/sum(spMap_cnt,"all");

%% linear weight adjustment
weight1_nuc = weight1(1);
weight1_cyto = weight1(2);
weight1_outcell = weight1(3);

weight_map1 = ones(size(spMap_cnt))*weight1_cyto;
I1 = sub2ind(size(weight_map1),coord_out_cell(:,2),coord_out_cell(:,1));
weight_map1(I1) = weight1_outcell;
I2 = sub2ind(size(weight_map1),coord_in_nuc(:,2),coord_in_nuc(:,1));
weight_map1(I2) = weight1_nuc;
clear I1 I2

%% power weight adjustment
weight2_nuc = weight2(1);
weight2_cyto = weight2(2);
weight2_outcell = weight2(3);

weight_map2 = ones(size(spMap_cnt))*weight2_outcell;
I1 = sub2ind(size(weight_map2),coord_in_cyto(:,2),coord_in_cyto(:,1));
weight_map2(I1) = weight2_cyto;
I2 = sub2ind(size(weight_map2),coord_in_nuc(:,2),coord_in_nuc(:,1));
weight_map2(I2) = weight2_nuc;
clear I1 I2

%% weight_matrix calculation by different methods
switch method
    case 'linear'
        weight_matrix = weight_by_cnt.*weight_map1;
    case 'power'
        weight_matrix = weight_by_cnt.^weight_map2;
    case 'combine'
        weight_matrix = (weight_by_cnt.*weight_map1).^weight_map2;
    otherwise
        disp('Method invalide, "linear" was used')
        weight_matrix = weight_by_cnt.*weight_map1;
end
weight_matrix = weight_matrix/sum(weight_matrix,'all');
end