%% summary_results.m
%  汇总所有数据集的DAC、AC和LoG-GP结果，并绘图
clc; close all;

datasets    = {'KIN40K', 'POL', 'PUMADYN32NM', 'SARCOS'};
dac_modes   = {'moe', 'gpoe', 'moe_ac', 'gpoe_ac'};
log_modes   = {'MOE', 'GPOE'};
all_labels  = {'LoG-MOE','LoG-GPOE','DAC-moe','DAC-gpoe','DAC-moe_ac','DAC-gpoe_ac'};
colors = [0.2 0.4 0.8;
          0.1 0.6 0.9;
          0.8 0.2 0.2;
          0.9 0.5 0.1;
          0.2 0.7 0.3;
          0.5 0.8 0.2];

metrics    = {'SMSE', 'RMSE', 'NLPD'};
N_datasets = numel(datasets);
N_methods  = numel(all_labels);
N_metrics  = numel(metrics);

% [SMSE, RMSE, NLPD, Train_T, Test_T]
all_results = nan(N_datasets, N_methods, 5);

%% 1. 加载结果
for di = 1:N_datasets
    ds = datasets{di};

    % LoG-GP
    f = fullfile('Result', 'Dataset', ds, 'LoGGP_results.mat');
    if exist(f, 'file')
        d = load(f);
        all_results(di, 1, :) = d.results(1,:);  % MOE
        all_results(di, 2, :) = d.results(2,:);  % GPOE
    end

    % DAC/AC
    for mi = 1:numel(dac_modes)
        f = fullfile('Result', 'Dataset', ds, [dac_modes{mi}, '.mat']);
        if exist(f, 'file')
            d = load(f);
            all_results(di, mi+2, :) = [d.smse, d.rmse, d.nlpd, d.t_train, d.t_test];
        end
    end
end

%% 2. 打印表格
for di = 1:N_datasets
    ds = datasets{di};
    fprintf('\n%s\n', repmat('=',1,76));
    fprintf('  Dataset: %s\n', ds);
    fprintf('%s\n', repmat('=',1,76));
    fprintf('  %-16s  %8s  %8s  %8s  %10s  %10s\n', 'Method','SMSE','RMSE','NLPD','Train_T(s)','Test_T(s)');
    fprintf('  %-16s  %8s  %8s  %8s  %10s  %10s\n', '------','----','----','----','----------','---------');

    for mi = 1:N_methods
        vals = squeeze(all_results(di, mi, :));
        if any(isnan(vals))
            fprintf('  %-16s  %8s  %8s  %8s  %10s  %10s\n', all_labels{mi}, 'N/A','N/A','N/A','N/A','N/A');
        else
            fprintf('  %-16s  %8.4f  %8.4f  %8.4f  %10.2f  %10.2f\n', ...
                all_labels{mi}, vals(1), vals(2), vals(3), vals(4), vals(5));
        end
    end
end
fprintf('\n%s\n', repmat('=',1,76));

%% 3. 绘图（每个指标一张，x轴数据集，每组bar是不同方法）
for ki = 1:N_metrics
    figure('Color','w','Position',[100+ki*50, 100, 1000, 420]);
    data_to_plot = squeeze(all_results(:, :, ki));  % N_datasets x N_methods
    b = bar(data_to_plot, 'grouped');
    for mi = 1:N_methods
        b(mi).FaceColor = colors(mi,:);
    end
    set(gca,'XTickLabel',datasets,'FontSize',11,'FontName','Times New Roman');
    ylabel(metrics{ki},'FontSize',13);
    title(sprintf('%s Comparison: LoG-GP vs DAC/AC', metrics{ki}), ...
        'FontSize',12,'FontName','Times New Roman');
    legend(all_labels,'Location','northwest','FontSize',9,'NumColumns',2);
    grid on; box on;
    saveas(gcf, fullfile('Result', sprintf('Summary_%s.png', metrics{ki})));
end

fprintf('完成。图片保存在 Result/ 目录。\n');