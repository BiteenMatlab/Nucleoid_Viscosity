function acc_fit = AccFit(accRange,Metric_optz,method)
    acc_range_fine = linspace(0,1,501);
    acc_fit_all = NaN(501,size(Metric_optz,2));
    for i = 1:size(Metric_optz,2)
        f = fit(accRange.',Metric_optz(:,i),'poly5');
        acc_fit_all(:,i) = f(acc_range_fine);
    end
    [~,acc_idx1] = min(acc_fit_all(:,1));
    [~,acc_idx2] = max(acc_fit_all(:,2));
    [~,acc_idx3] = max(acc_fit_all(:,3));
    acc_fit1 = acc_range_fine(acc_idx1);
    acc_fit2 = acc_range_fine(acc_idx2);
    acc_fit3 = acc_range_fine(acc_idx3);
    switch method
        case 'L2 Norm'
            acc_fit = acc_fit1;
        case 'PCC'
            acc_fit = acc_fit2;
        case 'SSIM'
            acc_fit = acc_fit3;
        case 'Combine'
            acc_fit = mean(rmoutliers([acc_fit1, acc_fit2, acc_fit3]));
        otherwise
            disp('Wrong method, Combine method is used')
            acc_fit = mean(rmoutliers([acc_fit1, acc_fit2, acc_fit3]));
    end
    fprintf('Nucleoid Accessibility = %.4f\n', acc_fit);

    if max(accRange) - min(accRange) > 0.5
        gap = 0.1;
    else
        gap = 0.02;
    end
    figure;
    figure('Units','pixels','Position',[100 100 800 200]);
    switch method
        case 'L2 Norm'
            plot(accRange, Metric_optz(:,1), 'LineWidth', 2);
            hold on
            plot(acc_range_fine, acc_fit_all(:,1), 'LineWidth', 1.5,'LineStyle','--','Color','r');
            hold on
            xl = xline(acc_fit,'-.'); 
            xl.LineWidth = 1.5;
            title('L2-Norm Distance (Minimize)', 'FontSize', 18);
            xticks(min(accRange):gap:max(accRange));
            set(gca, 'LineWidth', 1.5);
            xlabel('Accessibility')
        case 'PCC'
            plot(accRange, Metric_optz(:,2), 'LineWidth', 2);
            hold on
            plot(acc_range_fine, acc_fit_all(:,2), 'LineWidth', 1.5,'LineStyle','--','Color','r');
            hold on
            xl = xline(acc_fit,'-.'); 
            xl.LineWidth = 1.5;
            title('Pearson Correlation Coefficient (Maximize)', 'FontSize', 18);
            xticks(min(accRange):gap:max(accRange));
            set(gca, 'LineWidth', 1.5);
            xlabel('Accessibility')
        case 'SSIM'
            plot(accRange, Metric_optz(:,3), 'LineWidth', 2);
            hold on
            plot(acc_range_fine, acc_fit_all(:,3), 'LineWidth', 1.5,'LineStyle','--','Color','r');
            hold on
            xl = xline(acc_fit,'-.'); 
            xl.LineWidth = 1.5;
            title('SSIM (Maximize)', 'FontSize', 18);
            xticks(min(accRange):gap:max(accRange));
            set(gca, 'LineWidth', 1.5);
            xlabel('Accessibility')
        case 'Combine'
            figure('Units','pixels','Position',[100 100 1000 1500]);
            subplot(3, 1, 1)
            plot(accRange, Metric_optz(:,1), 'LineWidth', 2);
            hold on
            xl = xline(acc_fit,'-.'); 
            xl.LineWidth = 1.5;
            title('L2-Norm Distance (minimize)', 'FontSize', 18);
            xticks(min(accRange):gap:max(accRange));
            set(gca, 'LineWidth', 1.5);

            subplot(3, 1, 2)
            plot(accRange, Metric_optz(:,2), 'LineWidth', 2);
            hold on
            xl = xline(acc_fit,'-.'); 
            xl.LineWidth = 1.5;
            title('Pearson Correlation Coefficient (maximize)', 'FontSize', 18);
            xticks(min(accRange):gap:max(accRange));
            set(gca, 'LineWidth', 1.5);

            subplot(3, 1, 3)
            plot(accRange, Metric_optz(:,3), 'LineWidth', 2);
            hold on
            xl = xline(acc_fit,'-.'); 
            xl.LineWidth = 1.5;
            title('SSIM (maximize)', 'FontSize', 18);
            xticks(min(accRange):gap:max(accRange));
            set(gca, 'LineWidth', 1.5);
            xlabel('Accessibility')
    end
    
end