<<<<<<< HEAD
function nuc_outline_intp = Intp_Nuc(nuc_outline, num_pix_y, ratio)
% Intp_Nuc Interpolates the nucleoid outline to yield a smooth boundary.
%
% nuc_outline_intp = Intp_Nuc(nuc_outline, num_pix_y, ratio) interpolates the
% nucleoid outline (given as a set of (x,y) points in nuc_outline) by splitting
% the points into two groups (above and below the midline). For the upper side,
% it computes the minimum y (lower envelope) for each unique x value; for the
% lower side, it computes the maximum y (upper envelope). Each envelope is then
% interpolated using cubic spline interpolation at a resolution determined by 'ratio'.
%
% INPUTS:
% nuc_outline - A numeric [N x 2] matrix of nucleoid outline coordinates.
% (Column 1: x, Column 2: y.)
% num_pix_y - (Optional) Scalar specifying the number of divisions along the
% y-axis (default = 20).
% ratio - (Optional) Scalar specifying the interpolation resolution. A
% higher ratio yields more interpolated points (default = 20).
%
% OUTPUT:
% nuc_outline_intp - A [M x 2] matrix containing the interpolated nucleoid
% outline (each row is an (x,y) coordinate), where the
% first half corresponds to one envelope and the second half to the other.
%
% Example:
% % Given nuc_outline as Nx2 matrix:
% smooth_outline = Intp_Nuc(nuc_outline, 20, 20);
%
% Author: Xiaofeng Dai
% Date: 04/12/2025

%% 1. Set Defaults and Validate Inputs
if nargin < 2 || isempty(num_pix_y)
    num_pix_y = 20;
end
if nargin < 3 || isempty(ratio)
    ratio = 20;
end
if ~isnumeric(nuc_outline) || size(nuc_outline,2) ~= 2
    error('nuc_outline must be a numeric matrix with two columns.');
end

% Assume nuc_outline is given in normalized coordinate space,
% and we use num_pix_y to define a midline value.
midline = num_pix_y / 2;

%% 2. Split the Outline Points by Vertical Position
% Points with y > midline are assigned to group A (we will take the minimum y for each x)
% Points with y < midline are assigned to group B (we will take the maximum y for each x)
groupA = nuc_outline(nuc_outline(:,2) > midline, :);
groupB = nuc_outline(nuc_outline(:,2) < midline, :);

%% 3. Condense Each Group by X‐Value
% For groupA (typically the “lower” envelope of the upper half), for each unique x-value,
% take the minimum y-value.
groupA = sortrows(groupA, 1, 'ascend');
uniqueX_A = unique(groupA(:,1));
envelopeA = zeros(length(uniqueX_A),2);
for i = 1:length(uniqueX_A)
    envelopeA(i,1) = uniqueX_A(i);
    % For groupA, use min y.
    envelopeA(i,2) = min(groupA(groupA(:,1)==uniqueX_A(i), 2));
end

% For groupB (typically the “upper” envelope of the lower half), for each unique x-value,
% take the maximum y-value.
groupB = sortrows(groupB, 1, 'ascend');
uniqueX_B = unique(groupB(:,1));
envelopeB = zeros(length(uniqueX_B),2);
for i = 1:length(uniqueX_B)
    envelopeB(i,1) = uniqueX_B(i);
    % For groupB, use max y.
    envelopeB(i,2) = max(groupB(groupB(:,1)==uniqueX_B(i), 2));
end

%% 4. Spline Interpolation of Each Envelope
% Create a uniform query vector over the x-range of each envelope.
x_query_A = (min(envelopeA(:,1)) : 1/ratio : max(envelopeA(:,1)))';
x_query_B = (min(envelopeB(:,1)) : 1/ratio : max(envelopeB(:,1)))';

% Spline interpolation for each envelope:
y_interp_A = interp1(envelopeA(:,1), envelopeA(:,2), x_query_A, 'spline');
y_interp_B = interp1(envelopeB(:,1), envelopeB(:,2), x_query_B, 'spline');

% Combine interpolated envelopes as two rows: first row is envelopeA, second row is envelopeB.
% Then, concatenate horizontally and transpose so that each row corresponds to (x,y).
interpA = [x_query_A, y_interp_A];
interpB = [x_query_B, y_interp_B];

nuc_outline_intp = [interpA; interpB];

% % (Optional) Sort the final outline by the x coordinate if desired.
% nuc_outline_intp = sortrows(nuc_outline_intp, 1);
end
=======
function nuc_outline_intp = Intp_Nuc(nuc_outline, num_pix_y, ratio)
% Intp_Nuc Interpolates the nucleoid outline to yield a smooth boundary.
%
% nuc_outline_intp = Intp_Nuc(nuc_outline, num_pix_y, ratio) interpolates the
% nucleoid outline (given as a set of (x,y) points in nuc_outline) by splitting
% the points into two groups (above and below the midline). For the upper side,
% it computes the minimum y (lower envelope) for each unique x value; for the
% lower side, it computes the maximum y (upper envelope). Each envelope is then
% interpolated using cubic spline interpolation at a resolution determined by 'ratio'.
%
% INPUTS:
% nuc_outline - A numeric [N x 2] matrix of nucleoid outline coordinates.
% (Column 1: x, Column 2: y.)
% num_pix_y - (Optional) Scalar specifying the number of divisions along the
% y-axis (default = 20).
% ratio - (Optional) Scalar specifying the interpolation resolution. A
% higher ratio yields more interpolated points (default = 20).
%
% OUTPUT:
% nuc_outline_intp - A [M x 2] matrix containing the interpolated nucleoid
% outline (each row is an (x,y) coordinate), where the
% first half corresponds to one envelope and the second half to the other.
%
% Example:
% % Given nuc_outline as Nx2 matrix:
% smooth_outline = Intp_Nuc(nuc_outline, 20, 20);
%
% Author: Xiaofeng Dai
% Date: 04/12/2025

%% 1. Set Defaults and Validate Inputs
if nargin < 2 || isempty(num_pix_y)
    num_pix_y = 20;
end
if nargin < 3 || isempty(ratio)
    ratio = 20;
end
if ~isnumeric(nuc_outline) || size(nuc_outline,2) ~= 2
    error('nuc_outline must be a numeric matrix with two columns.');
end

% Assume nuc_outline is given in normalized coordinate space,
% and we use num_pix_y to define a midline value.
midline = num_pix_y / 2;

%% 2. Split the Outline Points by Vertical Position
% Points with y > midline are assigned to group A (we will take the minimum y for each x)
% Points with y < midline are assigned to group B (we will take the maximum y for each x)
groupA = nuc_outline(nuc_outline(:,2) > midline, :);
groupB = nuc_outline(nuc_outline(:,2) < midline, :);

%% 3. Condense Each Group by X‐Value
% For groupA (typically the “lower” envelope of the upper half), for each unique x-value,
% take the minimum y-value.
groupA = sortrows(groupA, 1, 'ascend');
uniqueX_A = unique(groupA(:,1));
envelopeA = zeros(length(uniqueX_A),2);
for i = 1:length(uniqueX_A)
    envelopeA(i,1) = uniqueX_A(i);
    % For groupA, use min y.
    envelopeA(i,2) = min(groupA(groupA(:,1)==uniqueX_A(i), 2));
end

% For groupB (typically the “upper” envelope of the lower half), for each unique x-value,
% take the maximum y-value.
groupB = sortrows(groupB, 1, 'ascend');
uniqueX_B = unique(groupB(:,1));
envelopeB = zeros(length(uniqueX_B),2);
for i = 1:length(uniqueX_B)
    envelopeB(i,1) = uniqueX_B(i);
    % For groupB, use max y.
    envelopeB(i,2) = max(groupB(groupB(:,1)==uniqueX_B(i), 2));
end

%% 4. Spline Interpolation of Each Envelope
% Create a uniform query vector over the x-range of each envelope.
x_query_A = (min(envelopeA(:,1)) : 1/ratio : max(envelopeA(:,1)))';
x_query_B = (min(envelopeB(:,1)) : 1/ratio : max(envelopeB(:,1)))';

% Spline interpolation for each envelope:
y_interp_A = interp1(envelopeA(:,1), envelopeA(:,2), x_query_A, 'spline');
y_interp_B = interp1(envelopeB(:,1), envelopeB(:,2), x_query_B, 'spline');

% Combine interpolated envelopes as two rows: first row is envelopeA, second row is envelopeB.
% Then, concatenate horizontally and transpose so that each row corresponds to (x,y).
interpA = [x_query_A, y_interp_A];
interpB = [x_query_B, y_interp_B];

nuc_outline_intp = [interpA; interpB];

% % (Optional) Sort the final outline by the x coordinate if desired.
% nuc_outline_intp = sortrows(nuc_outline_intp, 1);
end
>>>>>>> 5e4cf7b (Initial commit)
