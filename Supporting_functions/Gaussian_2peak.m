<<<<<<< HEAD
function params = Gaussian_2peak(distributionInput, numBins)
% fitDoubleGaussian Fits a probability distribution to a model consisting of two
% Gaussian curves plus a constant offset.
%
% params = fitDoubleGaussian(distributionInput) fits the normalized histogram
% of distributionInput (a numeric vector) to the model:
%
% f(x) = p1 * exp(-((x - p2)/p3)^2) + p4 * exp(-((x-p5)/p6)^2) + p7
%
% with default numBins = 50.
%
% params = fitDoubleGaussian(distributionInput, numBins) specifies the number
% of bins for the histogram.
%
% INPUT:
% distributionInput : Numeric vector of data to be analyzed.
% numBins : (Optional) Number of bins for histogram (default=50).
%
% OUTPUT:
% params : 1×7 vector with fitted parameters:
% p1: amplitude of Gaussian 1
% p2: center of Gaussian 1
% p3: sigma (width) of Gaussian 1
% p4: amplitude of Gaussian 2
% p5: center of Gaussian 2
% p6: sigma (width) of Gaussian 2
% p7: constant offset
%
% Additionally, the function generates a figure plotting the normalized
% histogram of the data, the two Gaussian components, and the residual between
% the data and the Gaussian model.
%
% Requires: Optimization Toolbox (lsqcurvefit).
%
% Example:
% params = fitDoubleGaussian(myData, 60);
%
% Author: Xiaofeng Dai
% Date: 03/28/2024

%% 1. Input Validation and Default Parameter Setup
narginchk(1,2);

if nargin < 2 || isempty(numBins)
    numBins = 50;
end

if ~isnumeric(distributionInput) || isempty(distributionInput)
    error('distributionInput must be a non-empty numeric vector.');
end

if ~isscalar(numBins) || ~isnumeric(numBins) || numBins <= 0
    error('numBins must be a positive numeric scalar.');
end

%% 2. Compute Histogram Data (Normalized) and Bin Midpoints
% Using histcounts to obtain normalized counts (i.e., probability).
[histCounts, histEdges] = histcounts(distributionInput, numBins, 'Normalization', 'probability');
yHistogram = histCounts(:);  % Ensure column vector

% Compute midpoints of each histogram bin.
xMidpoints = (histEdges(1:end-1) + histEdges(2:end)) / 2;
xMidpoints = xMidpoints.';

%% 3. Fit the Histogram Data to a Two-Gaussian Model
% The model is: 
%   f(x;p) = p1*exp(-((x-p2)/p3).^2) + p4*exp(-((x-p5)/p6).^2) + p7
%
% The helper function G2P uses lsqcurvefit to estimate parameters.
params = G2P(xMidpoints, yHistogram);

%% 4. Plot Data and Fitted Curves
figure;
hold on;

% Plot histogram of the input data (normalized as a probability density)
histogram(distributionInput, numBins, 'Normalization', 'probability', ...
          'FaceColor',[0.4660 0.6740 0.1880], 'EdgeColor','none');

% Calculate the two Gaussian components using fitted parameters.
gaussian1 = params(1) * exp(-((xMidpoints - params(2)) / params(3)).^2);
gaussian2 = params(4) * exp(-((xMidpoints - params(5)) / params(6)).^2);

% Compute the residual (difference between data and the sum of Gaussians).
residual = yHistogram - gaussian1 - gaussian2;

% Plot the fitted Gaussian components and residual.
plot(xMidpoints, gaussian1, 'r--', 'LineWidth', 1.5);
plot(xMidpoints, gaussian2, 'y--', 'LineWidth', 1.5);
plot(xMidpoints, residual,  'b--', 'LineWidth', 1.5);

% Enhance the plot with labels, legend, and title.
% xlabel('Data Value');
ylabel('Probability');
% title('Two-Gaussian Fit');
% legend({'Histogram','Gaussian 1','Gaussian 2','Residual'}, 'Location','Best');
hold off;
end

% Helper Function: G2P
%
% Purpose:
% Fit the data (xdata, ydata) to a dual Gaussian model plus constant offset.
%
% Model:
% f(x;p) = p(1)*exp(-((x-p(2))/p(3)).^2) + p(4)*exp(-((x-p(5))/p(6)).^2) + p(7)
%
% INPUT:
% xdata : Vector of independent variable values (e.g. bin midpoints).
% ydata : Column vector of dependent variable values (normalized histogram counts).
%
% OUTPUT:
% p_fit : 1×7 vector of fitted model parameters.
% err : Residual error from the least-squares optimization.
%
% Note: lsqcurvefit requires the Optimization Toolbox.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_fit, err] = G2P(xdata, ydata)
% Define the model function using an anonymous function.
modelFunc = @(p, x) p(1) * exp(-((x - p(2)) / p(3)).^2) + ...
                    p(4) * exp(-((x - p(5)) / p(6)).^2) + p(7);

% Set initial guesses for the parameters.
% p(1): amplitude of Gaussian 1
% p(2): center of Gaussian 1
% p(3): width (sigma) of Gaussian 1
% p(4): amplitude of Gaussian 2
% p(5): center of Gaussian 2
% p(6): width (sigma) of Gaussian 2
% p(7): constant offset
xMin = min(xdata);
xMax = max(xdata);
p_initial = [0.1, (xMax - xMin) / 4 + xMin, xMax / 20, ...
             0.1, 1.5 * (xMax - xMin) / 4 + xMin, xMax / 10, ...
             min(ydata) / 20];

% Lower bounds for the parameters
lb = [0, xMin, 0, 0, xMin, 0, 0];

% Upper bounds for the parameters
ub = [1, xMax, xMax, 1, xMax, xMax, max(ydata)];

% Set optimization options (suppress output)
options = optimset('Display', 'off');

% Perform least-squares curve fitting.
[p_fit, err] = lsqcurvefit(modelFunc, p_initial, xdata, ydata, lb, ub, options);
=======
function params = Gaussian_2peak(distributionInput, numBins)
% fitDoubleGaussian Fits a probability distribution to a model consisting of two
% Gaussian curves plus a constant offset.
%
% params = fitDoubleGaussian(distributionInput) fits the normalized histogram
% of distributionInput (a numeric vector) to the model:
%
% f(x) = p1 * exp(-((x - p2)/p3)^2) + p4 * exp(-((x-p5)/p6)^2) + p7
%
% with default numBins = 50.
%
% params = fitDoubleGaussian(distributionInput, numBins) specifies the number
% of bins for the histogram.
%
% INPUT:
% distributionInput : Numeric vector of data to be analyzed.
% numBins : (Optional) Number of bins for histogram (default=50).
%
% OUTPUT:
% params : 1×7 vector with fitted parameters:
% p1: amplitude of Gaussian 1
% p2: center of Gaussian 1
% p3: sigma (width) of Gaussian 1
% p4: amplitude of Gaussian 2
% p5: center of Gaussian 2
% p6: sigma (width) of Gaussian 2
% p7: constant offset
%
% Additionally, the function generates a figure plotting the normalized
% histogram of the data, the two Gaussian components, and the residual between
% the data and the Gaussian model.
%
% Requires: Optimization Toolbox (lsqcurvefit).
%
% Example:
% params = fitDoubleGaussian(myData, 60);
%
% Author: Xiaofeng Dai
% Date: 03/28/2024

%% 1. Input Validation and Default Parameter Setup
narginchk(1,2);

if nargin < 2 || isempty(numBins)
    numBins = 50;
end

if ~isnumeric(distributionInput) || isempty(distributionInput)
    error('distributionInput must be a non-empty numeric vector.');
end

if ~isscalar(numBins) || ~isnumeric(numBins) || numBins <= 0
    error('numBins must be a positive numeric scalar.');
end

%% 2. Compute Histogram Data (Normalized) and Bin Midpoints
% Using histcounts to obtain normalized counts (i.e., probability).
[histCounts, histEdges] = histcounts(distributionInput, numBins, 'Normalization', 'probability');
yHistogram = histCounts(:);  % Ensure column vector

% Compute midpoints of each histogram bin.
xMidpoints = (histEdges(1:end-1) + histEdges(2:end)) / 2;
xMidpoints = xMidpoints.';

%% 3. Fit the Histogram Data to a Two-Gaussian Model
% The model is: 
%   f(x;p) = p1*exp(-((x-p2)/p3).^2) + p4*exp(-((x-p5)/p6).^2) + p7
%
% The helper function G2P uses lsqcurvefit to estimate parameters.
params = G2P(xMidpoints, yHistogram);

%% 4. Plot Data and Fitted Curves
figure;
hold on;

% Plot histogram of the input data (normalized as a probability density)
histogram(distributionInput, numBins, 'Normalization', 'probability', ...
          'FaceColor',[0.4660 0.6740 0.1880], 'EdgeColor','none');

% Calculate the two Gaussian components using fitted parameters.
gaussian1 = params(1) * exp(-((xMidpoints - params(2)) / params(3)).^2);
gaussian2 = params(4) * exp(-((xMidpoints - params(5)) / params(6)).^2);

% Compute the residual (difference between data and the sum of Gaussians).
residual = yHistogram - gaussian1 - gaussian2;

% Plot the fitted Gaussian components and residual.
plot(xMidpoints, gaussian1, 'r--', 'LineWidth', 1.5);
plot(xMidpoints, gaussian2, 'y--', 'LineWidth', 1.5);
plot(xMidpoints, residual,  'b--', 'LineWidth', 1.5);

% Enhance the plot with labels, legend, and title.
% xlabel('Data Value');
ylabel('Probability');
% title('Two-Gaussian Fit');
% legend({'Histogram','Gaussian 1','Gaussian 2','Residual'}, 'Location','Best');
hold off;
end

% Helper Function: G2P
%
% Purpose:
% Fit the data (xdata, ydata) to a dual Gaussian model plus constant offset.
%
% Model:
% f(x;p) = p(1)*exp(-((x-p(2))/p(3)).^2) + p(4)*exp(-((x-p(5))/p(6)).^2) + p(7)
%
% INPUT:
% xdata : Vector of independent variable values (e.g. bin midpoints).
% ydata : Column vector of dependent variable values (normalized histogram counts).
%
% OUTPUT:
% p_fit : 1×7 vector of fitted model parameters.
% err : Residual error from the least-squares optimization.
%
% Note: lsqcurvefit requires the Optimization Toolbox.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_fit, err] = G2P(xdata, ydata)
% Define the model function using an anonymous function.
modelFunc = @(p, x) p(1) * exp(-((x - p(2)) / p(3)).^2) + ...
                    p(4) * exp(-((x - p(5)) / p(6)).^2) + p(7);

% Set initial guesses for the parameters.
% p(1): amplitude of Gaussian 1
% p(2): center of Gaussian 1
% p(3): width (sigma) of Gaussian 1
% p(4): amplitude of Gaussian 2
% p(5): center of Gaussian 2
% p(6): width (sigma) of Gaussian 2
% p(7): constant offset
xMin = min(xdata);
xMax = max(xdata);
p_initial = [0.1, (xMax - xMin) / 4 + xMin, xMax / 20, ...
             0.1, 1.5 * (xMax - xMin) / 4 + xMin, xMax / 10, ...
             min(ydata) / 20];

% Lower bounds for the parameters
lb = [0, xMin, 0, 0, xMin, 0, 0];

% Upper bounds for the parameters
ub = [1, xMax, xMax, 1, xMax, xMax, max(ydata)];

% Set optimization options (suppress output)
options = optimset('Display', 'off');

% Perform least-squares curve fitting.
[p_fit, err] = lsqcurvefit(modelFunc, p_initial, xdata, ydata, lb, ub, options);
>>>>>>> 5e4cf7b (Initial commit)
end