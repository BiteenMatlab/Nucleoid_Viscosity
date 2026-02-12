<<<<<<< HEAD
function [BtstrpResult, ind_pass]=Btstrp_showResult(Btstrp_NucVol,Btstrp_Acc,Btstrp_Vis_Nuc,Btstrp_Vis_Cyto,Btstrp_Vis_PeriCore,Btstrp_Vis_Domain,method)
    if nargin<7 || isempty(method)
        method = {'median','on'};
    end
    if ~isempty(method) && ischar(method)
        method = {method, 'on'};
    end


    Book_all = cat(2,Btstrp_NucVol,Btstrp_Acc,Btstrp_Vis_Nuc,Btstrp_Vis_Cyto);
    if ~isempty(Btstrp_Vis_PeriCore)
        Book_all = cat(2, Book_all, Btstrp_Vis_PeriCore);
    end
    if ~isempty(Btstrp_Vis_Domain)
        Book_all = cat(2, Book_all, Btstrp_Vis_Domain);
    end
    Deter = zeros(size(Book_all));
    %%
    if ~strcmpi(method{1},'off')
        for i = 1:size(Book_all,2)
            deter = isoutlier(Book_all(:,i),method{1})|isnan(Book_all(:,i));
            Deter(:,i) = deter;
        end
    clear i deter
    end
    %%
    if strcmpi(method{2},'on') && ~isempty(Btstrp_Vis_PeriCore)
        deter = zeros(size(Btstrp_Vis_Nuc));
        for i = 1:size(Btstrp_Vis_Nuc,1)
            if Btstrp_Vis_Nuc(i)<min(Btstrp_Vis_PeriCore(i,:)) || Btstrp_Vis_Nuc(i)>max(Btstrp_Vis_PeriCore(i,:))
                deter(i) = 1;
            end
        end
        Deter = cat(2,Deter,deter);
        clear deter i 
    end
    %%
    ind_pass = find(sum(Deter,2)==0);
    clear Deter Book_all
    %%
    iter_gap = 1;
    % Nucleoid volume ratio and accessibility
    [NucVol_mean,NucVol_std]=...
        Result_show(...
        Btstrp_NucVol(ind_pass,:),...
        iter_gap,'on','on','Volume Ratio');
    [Acc_mean,Acc_std]=...
        Result_show(...
        Btstrp_Acc(ind_pass,:),...
        iter_gap,'on','on','Accessibility');
    disp('================================================')
    % Nucleoid and cytoplasm viscosity
    [VisNuc_mean,VisNuc_std]=...
        Result_show(...
        Btstrp_Vis_Nuc(ind_pass,:),...
        iter_gap,'on','on','Nucleoid Viscosity');
    [VisCyto_mean,VisCyto_std]=...
        Result_show(...
        Btstrp_Vis_Cyto(ind_pass,:),...
        iter_gap,'on','on','Cytoplasm Viscosity');
    disp('================================================')

    if ~isempty(Btstrp_Vis_PeriCore)
    % Viscosity of nucleoid periphery and core
    [VisPeri_mean,VisPeri_std]=...
        Result_show(...
        Btstrp_Vis_PeriCore(ind_pass,1),...
        iter_gap,'on','on','Periphery Viscosity');
    [VisCore_mean,VisCore_std]=...
        Result_show(...
        Btstrp_Vis_PeriCore(ind_pass,2),...
        iter_gap,'on','on','Core Viscosity');
    disp('================================================')
    end

    if ~isempty(Btstrp_Vis_Domain)
    % Viscosity of macrodomain
    Domain_result = nan(size(Btstrp_Vis_Domain,2),2);
    for i = 1:size(Btstrp_Vis_Domain,2)
        [VisDomain_mean,VisDomain_std]=...
            Result_show(...
            Btstrp_Vis_Domain(ind_pass,i),...
            iter_gap,'on','on',strcat('Macrodomain_',string(i)));
        Domain_result(i,:) = [VisDomain_mean,VisDomain_std];
        disp('-------------------------------------------------')
    end
    end

    BtstrpResult.NucVol  = [NucVol_mean,NucVol_std];
    BtstrpResult.NucAcc  = [Acc_mean,Acc_std];
    BtstrpResult.NucVis  = [VisNuc_mean,VisNuc_std];
    BtstrpResult.CytoVis = [VisCyto_mean,VisCyto_std];
    if ~isempty(Btstrp_Vis_PeriCore)
        BtstrpResult.NucPeriCore = cat(1,[VisPeri_mean,VisPeri_std],[VisCore_mean,VisCore_std]);
    end
    if ~isempty(Btstrp_Vis_Domain)
        BtstrpResult.NucDomain = Domain_result;
    end
end



function [Mean, Std, figHist, figBootstrap] =...
    Result_show(bootstrap, iterGap, figPlotHist, figPlotBootstrap,Title)
%% ----- Input Handling and Defaults -----
if nargin < 5 || isempty(Title)
    Title = '';
end
if nargin < 4 || isempty(figPlotBootstrap)
    figPlotBootstrap = 'off';
end
if nargin < 3 || isempty(figPlotHist)
    figPlotHist = 'on';
end
if nargin < 2 || isempty(iterGap)
    iterGap = 5;
end

% Defensive check: bootstrap data must be numeric vectors.
if ~isnumeric(bootstrap) || ~isvector(bootstrap)
    error('Bootstrapping result1 must be a numeric vector.');
end

%% ----- Compute Statistics from the Cleaned Data -----
Mean = mean(bootstrap);
Std  = std(bootstrap);


%% ----- Create Histogram Figure -----
figHist = figure;
figure('Units','pixels','Position',[100 100 800 500]);
% histogram(cleanNucVol,...
%     linspace(min(bootstrap)-0.2*(max(bootstrap)-min(bootstrap)),...
%     max(bootstrap)+0.2*(max(bootstrap)-min(bootstrap)),10))
h = histogram(bootstrap);
ylim([0, max(h.Values) + 1]);
histTitleNuc = sprintf('%.2f ± %.2f', mean(bootstrap), std(bootstrap));
title(Title, histTitleNuc, 'FontSize', 15)
ylabel('Counts');
% Close the histogram figure if the plotting flag is set to 'off'.
if strcmpi(figPlotHist, 'off')
    close(figHist)
end


%% ----- Compute Bootstrap Convergence Trends -----
% Create an array of indices stepping by iterGap until the number of cleaned samples.
iterationIdx = iterGap:iterGap:length(bootstrap);
numIter = numel(iterationIdx);

% Preallocate arrays for cumulative means and standard deviations.
bsMean = zeros(numIter, 1);
bsStd  = zeros(numIter, 1);

% Compute cumulative statistics.
for i = 1:numIter
    currentData = bootstrap(1:iterationIdx(i));
    bsMean(i) = mean(currentData);
    bsStd(i)  = std(currentData);
end

% Create a new figure to display bootstrap convergence for property 1.
figBootstrap = figure;
figure('Units','pixels','Position',[100 100 800 500]);
% Plot cumulative mean for property 1.
subplot(1,2,1)
plot(iterationIdx, bsMean, 'LineWidth', 1.5)
title('Mean Value', 'FontWeight', 'normal')
ylabel(Title)
xlabel('Bootstrap Iteration')

% Plot cumulative standard deviation for property 1.
subplot(1,2,2)
plot(iterationIdx, bsStd, 'LineWidth', 1.5, 'Color', '#EB9330')
title('Standard Deviation', 'FontWeight', 'normal')
% ylabel(Title{1})
xlabel('Bootstrap Iteration')
sgtitle('Mean and Std vs. Bootstrap Iterations', 'FontWeight', 'bold')
% figBootstrap.Position = [100 100 1000 750];
% Close the bootstrap figure if the corresponding flag is set to 'off'.
if strcmpi(figPlotBootstrap, 'off')
    close(figBootstrap)
end
end
=======
function [BtstrpResult, ind_pass]=Btstrp_showResult(Btstrp_NucVol,Btstrp_Acc,Btstrp_Vis_Nuc,Btstrp_Vis_Cyto,Btstrp_Vis_PeriCore,Btstrp_Vis_Domain,method)
    if nargin<7 || isempty(method)
        method = {'median','on'};
    end
    if ~isempty(method) && ischar(method)
        method = {method, 'on'};
    end


    Book_all = cat(2,Btstrp_NucVol,Btstrp_Acc,Btstrp_Vis_Nuc,Btstrp_Vis_Cyto);
    if ~isempty(Btstrp_Vis_PeriCore)
        Book_all = cat(2, Book_all, Btstrp_Vis_PeriCore);
    end
    if ~isempty(Btstrp_Vis_Domain)
        Book_all = cat(2, Book_all, Btstrp_Vis_Domain);
    end
    Deter = zeros(size(Book_all));
    %%
    if ~strcmpi(method{1},'off')
        for i = 1:size(Book_all,2)
            deter = isoutlier(Book_all(:,i),method{1})|isnan(Book_all(:,i));
            Deter(:,i) = deter;
        end
    clear i deter
    end
    %%
    if strcmpi(method{2},'on') && ~isempty(Btstrp_Vis_PeriCore)
        deter = zeros(size(Btstrp_Vis_Nuc));
        for i = 1:size(Btstrp_Vis_Nuc,1)
            if Btstrp_Vis_Nuc(i)<min(Btstrp_Vis_PeriCore(i,:)) || Btstrp_Vis_Nuc(i)>max(Btstrp_Vis_PeriCore(i,:))
                deter(i) = 1;
            end
        end
        Deter = cat(2,Deter,deter);
        clear deter i 
    end
    %%
    ind_pass = find(sum(Deter,2)==0);
    clear Deter Book_all
    %%
    iter_gap = 1;
    % Nucleoid volume ratio and accessibility
    [NucVol_mean,NucVol_std]=...
        Result_show(...
        Btstrp_NucVol(ind_pass,:),...
        iter_gap,'on','on','Volume Ratio');
    [Acc_mean,Acc_std]=...
        Result_show(...
        Btstrp_Acc(ind_pass,:),...
        iter_gap,'on','on','Accessibility');
    disp('================================================')
    % Nucleoid and cytoplasm viscosity
    [VisNuc_mean,VisNuc_std]=...
        Result_show(...
        Btstrp_Vis_Nuc(ind_pass,:),...
        iter_gap,'on','on','Nucleoid Viscosity');
    [VisCyto_mean,VisCyto_std]=...
        Result_show(...
        Btstrp_Vis_Cyto(ind_pass,:),...
        iter_gap,'on','on','Cytoplasm Viscosity');
    disp('================================================')

    if ~isempty(Btstrp_Vis_PeriCore)
    % Viscosity of nucleoid periphery and core
    [VisPeri_mean,VisPeri_std]=...
        Result_show(...
        Btstrp_Vis_PeriCore(ind_pass,1),...
        iter_gap,'on','on','Periphery Viscosity');
    [VisCore_mean,VisCore_std]=...
        Result_show(...
        Btstrp_Vis_PeriCore(ind_pass,2),...
        iter_gap,'on','on','Core Viscosity');
    disp('================================================')
    end

    if ~isempty(Btstrp_Vis_Domain)
    % Viscosity of macrodomain
    Domain_result = nan(size(Btstrp_Vis_Domain,2),2);
    for i = 1:size(Btstrp_Vis_Domain,2)
        [VisDomain_mean,VisDomain_std]=...
            Result_show(...
            Btstrp_Vis_Domain(ind_pass,i),...
            iter_gap,'on','on',strcat('Macrodomain_',string(i)));
        Domain_result(i,:) = [VisDomain_mean,VisDomain_std];
        disp('-------------------------------------------------')
    end
    end

    BtstrpResult.NucVol  = [NucVol_mean,NucVol_std];
    BtstrpResult.NucAcc  = [Acc_mean,Acc_std];
    BtstrpResult.NucVis  = [VisNuc_mean,VisNuc_std];
    BtstrpResult.CytoVis = [VisCyto_mean,VisCyto_std];
    if ~isempty(Btstrp_Vis_PeriCore)
        BtstrpResult.NucPeriCore = cat(1,[VisPeri_mean,VisPeri_std],[VisCore_mean,VisCore_std]);
    end
    if ~isempty(Btstrp_Vis_Domain)
        BtstrpResult.NucDomain = Domain_result;
    end
end



function [Mean, Std, figHist, figBootstrap] =...
    Result_show(bootstrap, iterGap, figPlotHist, figPlotBootstrap,Title)
%% ----- Input Handling and Defaults -----
if nargin < 5 || isempty(Title)
    Title = '';
end
if nargin < 4 || isempty(figPlotBootstrap)
    figPlotBootstrap = 'off';
end
if nargin < 3 || isempty(figPlotHist)
    figPlotHist = 'on';
end
if nargin < 2 || isempty(iterGap)
    iterGap = 5;
end

% Defensive check: bootstrap data must be numeric vectors.
if ~isnumeric(bootstrap) || ~isvector(bootstrap)
    error('Bootstrapping result1 must be a numeric vector.');
end

%% ----- Compute Statistics from the Cleaned Data -----
Mean = mean(bootstrap);
Std  = std(bootstrap);


%% ----- Create Histogram Figure -----
figHist = figure;
figure('Units','pixels','Position',[100 100 800 500]);
% histogram(cleanNucVol,...
%     linspace(min(bootstrap)-0.2*(max(bootstrap)-min(bootstrap)),...
%     max(bootstrap)+0.2*(max(bootstrap)-min(bootstrap)),10))
h = histogram(bootstrap);
ylim([0, max(h.Values) + 1]);
histTitleNuc = sprintf('%.2f ± %.2f', mean(bootstrap), std(bootstrap));
title(Title, histTitleNuc, 'FontSize', 15)
ylabel('Counts');
% Close the histogram figure if the plotting flag is set to 'off'.
if strcmpi(figPlotHist, 'off')
    close(figHist)
end


%% ----- Compute Bootstrap Convergence Trends -----
% Create an array of indices stepping by iterGap until the number of cleaned samples.
iterationIdx = iterGap:iterGap:length(bootstrap);
numIter = numel(iterationIdx);

% Preallocate arrays for cumulative means and standard deviations.
bsMean = zeros(numIter, 1);
bsStd  = zeros(numIter, 1);

% Compute cumulative statistics.
for i = 1:numIter
    currentData = bootstrap(1:iterationIdx(i));
    bsMean(i) = mean(currentData);
    bsStd(i)  = std(currentData);
end

% Create a new figure to display bootstrap convergence for property 1.
figBootstrap = figure;
figure('Units','pixels','Position',[100 100 800 500]);
% Plot cumulative mean for property 1.
subplot(1,2,1)
plot(iterationIdx, bsMean, 'LineWidth', 1.5)
title('Mean Value', 'FontWeight', 'normal')
ylabel(Title)
xlabel('Bootstrap Iteration')

% Plot cumulative standard deviation for property 1.
subplot(1,2,2)
plot(iterationIdx, bsStd, 'LineWidth', 1.5, 'Color', '#EB9330')
title('Standard Deviation', 'FontWeight', 'normal')
% ylabel(Title{1})
xlabel('Bootstrap Iteration')
sgtitle('Mean and Std vs. Bootstrap Iterations', 'FontWeight', 'bold')
% figBootstrap.Position = [100 100 1000 750];
% Close the bootstrap figure if the corresponding flag is set to 'off'.
if strcmpi(figPlotBootstrap, 'off')
    close(figBootstrap)
end
end
>>>>>>> 5e4cf7b (Initial commit)
