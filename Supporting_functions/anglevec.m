function theta = anglevec(x, y)
% anglevec Computes the angle between two 2D vectors
%
% theta = anglevec(x, y) returns the angle (in radians) between the two
% 2-D vectors x and y. The angle is computed using the cosine of the angle:
%
% cos_val = dot(x, y) / (norm(x) * norm(y))
% base_angle = acos(cos_val)
%
% Since acos returns values only in [0, pi], the function then adjusts the
% result based on the sign of y(2) (and y(1) if y(2)==0) so that the output
% angle spans the full range [0, 2pi). This method assumes that vector x is used as
% a reference direction (e.g., if x is along the positive x-axis), and that the
% quadrant of y can be determined by its second component.
%
% Inputs:
% x - A 2-element numeric vector representing the first vector.
% y - A 2-element numeric vector representing the second vector.
%
% Output:
% theta - A scalar in [0, 2pi) representing the angle (in radians) from vector x to y.
%
% Defensive Mechanisms:
% - The function verifies that x and y are two-element vectors.
% - It checks that neither vector is the zero vector.
% - It clamps the computed cos_val to the interval [-1,1] to avoid
% potential numerical issues with acos.
%
% Example:
% % Compute angle between [1,0] and [0,1]:
% theta = anglevec([1 0], [0 1]); % returns pi/2
%
% Author: Xiaofeng Dai
% Date: 07/05/2023

%% 1. Input Validation
if ~isvector(x) || ~isvector(y) || numel(x) ~= 2 || numel(y) ~= 2
    error('Both x and y must be 2-element vectors.');
end
if norm(x) == 0 || norm(y) == 0
    error('Input vectors x and y must be nonzero.');
end

%% 2. Compute the Cosine of the Angle and Clamp it
cos_val = dot(x, y) / (norm(x) * norm(y));
% Clamp cos_val to the interval [-1, 1] to avoid acos domain errors due to round-off:
cos_val = max(min(cos_val, 1), -1);

% Compute the base angle (in [0,pi]) from the acos value.
base_angle = acos(cos_val);

%% 3. Determine the Correct Quadrant to Recover a Full [0, 2*pi) Angle
if y(2) < 0
    theta = 2*pi - base_angle;
elseif y(2) > 0
    theta = base_angle;
else
    if y(1) >= 0
        theta = 0;
    else % if y(1) < 0
        theta = pi;
    end
end
end
