<<<<<<< HEAD
function orgIndex = reorgIndex(Index)
% reorgIndex Reorganizes a two‐column numeric matrix into a cell array.
%
% orgIndex = reorgIndex(Index) takes a numeric matrix "Index" where the first
% column represents ROI identifiers and the second column represents cell IDs.
% It returns a cell array "orgIndex" whose length is equal to the maximum ROI
% identifier. For each ROI i, orgIndex{i} contains all cell IDs (from the second
% column of Index) corresponding to ROI i.
%
% Example:
% Index = [ 1 101;
% 2 201;
% 1 102;
% 2 202;
% 3 301];
% orgIndex = reorgIndex(Index);
% % orgIndex{1} returns [101;102], orgIndex{2} returns [201;202], etc.
%
% Note: If some ROI indices between 1 and max(Index(:,1)) are absent, the
% corresponding cell in orgIndex will be empty.
%
% Author: Xiaofeng Dai
% Date: 04/11/2025

%% Input Validation
if ~ismatrix(Index) || size(Index,2) < 2
    error('Input "Index" must be a numeric matrix with at least two columns.');
end

%% Sort the Input for Consistency
sortedIndex = sortrows(Index, 1);
maxROI = max(sortedIndex(:,1));

%% Preallocate Output Cell Array
orgIndex = cell(maxROI, 1);

%% Fill Each Cell with the Corresponding Cell IDs
for roi = 1:maxROI
    index_each_roi = sortedIndex(sortedIndex(:,1) == roi, 2);
    orgIndex{roi} = sort(index_each_roi);
end
=======
function orgIndex = reorgIndex(Index)
% reorgIndex Reorganizes a two‐column numeric matrix into a cell array.
%
% orgIndex = reorgIndex(Index) takes a numeric matrix "Index" where the first
% column represents ROI identifiers and the second column represents cell IDs.
% It returns a cell array "orgIndex" whose length is equal to the maximum ROI
% identifier. For each ROI i, orgIndex{i} contains all cell IDs (from the second
% column of Index) corresponding to ROI i.
%
% Example:
% Index = [ 1 101;
% 2 201;
% 1 102;
% 2 202;
% 3 301];
% orgIndex = reorgIndex(Index);
% % orgIndex{1} returns [101;102], orgIndex{2} returns [201;202], etc.
%
% Note: If some ROI indices between 1 and max(Index(:,1)) are absent, the
% corresponding cell in orgIndex will be empty.
%
% Author: Xiaofeng Dai
% Date: 04/11/2025

%% Input Validation
if ~ismatrix(Index) || size(Index,2) < 2
    error('Input "Index" must be a numeric matrix with at least two columns.');
end

%% Sort the Input for Consistency
sortedIndex = sortrows(Index, 1);
maxROI = max(sortedIndex(:,1));

%% Preallocate Output Cell Array
orgIndex = cell(maxROI, 1);

%% Fill Each Cell with the Corresponding Cell IDs
for roi = 1:maxROI
    index_each_roi = sortedIndex(sortedIndex(:,1) == roi, 2);
    orgIndex{roi} = sort(index_each_roi);
end
>>>>>>> 5e4cf7b (Initial commit)
end