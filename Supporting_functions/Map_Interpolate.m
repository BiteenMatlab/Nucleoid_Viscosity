<<<<<<< HEAD
function map_intp = Map_Interpolate(map, scale_rate, intp_method)
% Map_Interpolate Interpolates a 2D map to a finer grid.
%
% map_intp = Map_Interpolate(map) uses the default scale_rate of 10 and
% cubic interpolation.
%
% map_intp = Map_Interpolate(map, scale_rate) uses the specified scale_rate.
%
% map_intp = Map_Interpolate(map, scale_rate, intp_method) uses the specified
% interpolation method, such as 'cubic', 'linear', etc.
%
% The function assumes that the input map spans the unit square [0,1]×[0,1]
% with bins defined by equally spaced edges. It first computes the bin centers
% for the original grid, then builds a query grid that is scale_rate times finer
% along each dimension. It interpolates the map onto the new grid using interp2
% and scales the result by 1/(scale_rate^2) so that the map remains normalized.
%
% INPUTS:
% map - Numeric [m x n] matrix.
% scale_rate - (Optional) Scalar, the interpolation scale factor (default 10).
% intp_method - (Optional) String specifying the interp2 method (default 'cubic').
%
% OUTPUT:
% map_intp - The interpolated and rescaled map.
%
% Example:
% originalMap = rand(20,30);
% fineMap = Map_Interpolate(originalMap, 10, 'cubic');
%
% Author: Xiaofeng Dai
% Date: 04/12/2025

%% 1. Set default parameters
if nargin < 2 || isempty(scale_rate)
    scale_rate = 10;
end
if nargin < 3 || isempty(intp_method)
    intp_method = 'cubic';
end

%% 2. Input validation
if ~ismatrix(map) || ~isnumeric(map)
    error('Input "map" must be a numeric 2D matrix.');
end

%% 3. Define original grid based on map dimensions assuming unit interval
[num_rows, num_cols] = size(map);
dx = 1 / num_cols;
dy = 1 / num_rows;

% Define bin edges and then compute bin centers.
x_edges = 0:dx:1;
y_edges = 0:dy:1;
x = Get_bincenter(x_edges);
y = Get_bincenter(y_edges);
[X, Y] = meshgrid(x, y);

%% 4. Define query grid (scale_rate times finer)
dx_fine = dx / scale_rate;
dy_fine = dy / scale_rate;
xq_edges = 0:dx_fine:1;
yq_edges = 0:dy_fine:1;
xq = Get_bincenter(xq_edges);
yq = Get_bincenter(yq_edges);
[Xq, Yq] = meshgrid(xq, yq);

%% 5. Interpolate the map onto the fine grid
map_intp = interp2(X, Y, map, Xq, Yq, intp_method);

% Scale the interpolated map to account for the increased resolution.
map_intp = map_intp ./ (scale_rate^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Bincenter = Get_bincenter(edgeArray)
% Get_bincenter Computes the center of bins from an array of edges.
%
% Bincenter = Get_bincenter(edgeArray) returns a column vector containing the
% arithmetic mean of each pair of consecutive values in edgeArray.
%
% Example:
% edges = 0:0.1:1;
% centers = Get_bincenter(edges);
%
% INPUT:
% edgeArray - A numeric vector with at least two elements.
%
% OUTPUT:
% Bincenter - A column vector of bin centers.

if ~isnumeric(edgeArray) || numel(edgeArray) < 2
    error('edgeArray must be a numeric vector with at least two elements.');
end

% Vectorized computation of bin centers.
Bincenter = (edgeArray(1:end-1) + edgeArray(2:end)) / 2;
Bincenter = Bincenter(:);
=======
function map_intp = Map_Interpolate(map, scale_rate, intp_method)
% Map_Interpolate Interpolates a 2D map to a finer grid.
%
% map_intp = Map_Interpolate(map) uses the default scale_rate of 10 and
% cubic interpolation.
%
% map_intp = Map_Interpolate(map, scale_rate) uses the specified scale_rate.
%
% map_intp = Map_Interpolate(map, scale_rate, intp_method) uses the specified
% interpolation method, such as 'cubic', 'linear', etc.
%
% The function assumes that the input map spans the unit square [0,1]×[0,1]
% with bins defined by equally spaced edges. It first computes the bin centers
% for the original grid, then builds a query grid that is scale_rate times finer
% along each dimension. It interpolates the map onto the new grid using interp2
% and scales the result by 1/(scale_rate^2) so that the map remains normalized.
%
% INPUTS:
% map - Numeric [m x n] matrix.
% scale_rate - (Optional) Scalar, the interpolation scale factor (default 10).
% intp_method - (Optional) String specifying the interp2 method (default 'cubic').
%
% OUTPUT:
% map_intp - The interpolated and rescaled map.
%
% Example:
% originalMap = rand(20,30);
% fineMap = Map_Interpolate(originalMap, 10, 'cubic');
%
% Author: Xiaofeng Dai
% Date: 04/12/2025

%% 1. Set default parameters
if nargin < 2 || isempty(scale_rate)
    scale_rate = 10;
end
if nargin < 3 || isempty(intp_method)
    intp_method = 'cubic';
end

%% 2. Input validation
if ~ismatrix(map) || ~isnumeric(map)
    error('Input "map" must be a numeric 2D matrix.');
end

%% 3. Define original grid based on map dimensions assuming unit interval
[num_rows, num_cols] = size(map);
dx = 1 / num_cols;
dy = 1 / num_rows;

% Define bin edges and then compute bin centers.
x_edges = 0:dx:1;
y_edges = 0:dy:1;
x = Get_bincenter(x_edges);
y = Get_bincenter(y_edges);
[X, Y] = meshgrid(x, y);

%% 4. Define query grid (scale_rate times finer)
dx_fine = dx / scale_rate;
dy_fine = dy / scale_rate;
xq_edges = 0:dx_fine:1;
yq_edges = 0:dy_fine:1;
xq = Get_bincenter(xq_edges);
yq = Get_bincenter(yq_edges);
[Xq, Yq] = meshgrid(xq, yq);

%% 5. Interpolate the map onto the fine grid
map_intp = interp2(X, Y, map, Xq, Yq, intp_method);

% Scale the interpolated map to account for the increased resolution.
map_intp = map_intp ./ (scale_rate^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Bincenter = Get_bincenter(edgeArray)
% Get_bincenter Computes the center of bins from an array of edges.
%
% Bincenter = Get_bincenter(edgeArray) returns a column vector containing the
% arithmetic mean of each pair of consecutive values in edgeArray.
%
% Example:
% edges = 0:0.1:1;
% centers = Get_bincenter(edges);
%
% INPUT:
% edgeArray - A numeric vector with at least two elements.
%
% OUTPUT:
% Bincenter - A column vector of bin centers.

if ~isnumeric(edgeArray) || numel(edgeArray) < 2
    error('edgeArray must be a numeric vector with at least two elements.');
end

% Vectorized computation of bin centers.
Bincenter = (edgeArray(1:end-1) + edgeArray(2:end)) / 2;
Bincenter = Bincenter(:);
>>>>>>> 5e4cf7b (Initial commit)
end