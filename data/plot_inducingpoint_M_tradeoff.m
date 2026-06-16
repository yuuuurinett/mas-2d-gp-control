% plot_all_tradeoffs_batch_and_markdown.m
clear; clc; close all;

% =========================================================================
% ⚙️ 全局配置区
% =========================================================================
%datasets = {'KIN40K', 'POL', 'PUMADYN32NM', 'SARCOS'}; 
datasets = {'KIN40K'}; 
Method = 'poe';
train_ratio = 0.4;
tr_tag = round(train_ratio * 100);

M_list = [500, 1000, 1500, 2000, 2500];
seeds = 1:3;

all_table_rows = {}; 

fprintf('=================================================\n');
fprintf('🎨 开始批量绘图并生成终极 Markdown 排版...\n');
fprintf('=================================================\n');

for d_idx = 1:numel(datasets)
    DatasetName = datasets{d_idx};
    fprintf('>> 正在处理图表: %s...\n', DatasetName);
    
    ResultFolder = fullfile('result', 'Dataset', DatasetName);
    
    SMSE_all   = nan(numel(M_list), numel(seeds));
    TrainT_all = nan(numel(M_list), numel(seeds));
    TestT_all  = nan(numel(M_list), numel(seeds));
    
    for mi = 1:numel(M_list)
        M = M_list(mi);
        for si = 1:numel(seeds)
            seed = seeds(si);
            file_name = sprintf('%s_M%d_tr%d_mc%d.mat', Method, M, tr_tag, seed);
            file_path = fullfile(ResultFolder, file_name);

            if exist(file_path, 'file')
                S = load(file_path);
                if isfield(S, 'smse'), SMSE_all(mi, si) = S.smse; end
                
                if isfield(S, 't_train_per_point')
                    TrainT_all(mi, si) = S.t_train_per_point;   
                elseif isfield(S, 't_train_total') && isfield(S, 'N_train')
                    TrainT_all(mi, si) = S.t_train_total / S.N_train * 1000;
                end
                
                if isfield(S, 't_test_per_point')
                    TestT_all(mi, si) = S.t_test_per_point;     
                elseif isfield(S, 't_test_total') && isfield(S, 'N_eval')
                    TestT_all(mi, si) = S.t_test_total / S.N_eval * 1000;
                end
            end
        end
    end

    SMSE_mean   = mean(SMSE_all,   2, 'omitnan');
    SMSE_std    = std(SMSE_all,  0, 2, 'omitnan');
    TrainT_mean = mean(TrainT_all, 2, 'omitnan');
    TrainT_std  = std(TrainT_all,0, 2, 'omitnan');
    TestT_mean  = mean(TestT_all,  2, 'omitnan');
    TestT_std   = std(TestT_all, 0, 2, 'omitnan');

    switch upper(DatasetName)
        case 'KIN40K',      target_M = 2500;  
        case 'POL',         target_M = 2000;  
        case 'PUMADYN32NM', target_M = 500;
        case 'SARCOS',      target_M = 1500;  
        otherwise,          target_M = [];    
    end

    for i = 1:length(M_list)
        if M_list(i) == target_M
            b = '**'; 
        else
            b = '';   
        end
        
        if i == 1
            ds_name = sprintf('**%s**', DatasetName);
        else
            ds_name = '';
        end
        
        tr_str   = sprintf('%.4f±%.4f', TrainT_mean(i), TrainT_std(i));
        te_str   = sprintf('%.4f±%.4f', TestT_mean(i), TestT_std(i));
        smse_str = sprintf('%.6f±%.6f', SMSE_mean(i), SMSE_std(i));
        
        row_str = sprintf('| %s | %s%d%s | %s%s%s | %s%s%s | %s%s%s |', ...
            ds_name, b, M_list(i), b, b, tr_str, b, b, te_str, b, b, smse_str, b);
        all_table_rows{end+1} = row_str;
    end
    
    % (已删除空行)

    % 画图部分 (保持不变，省略注释)
    fig = figure('Name', sprintf('%s_Tradeoff', DatasetName), 'Color', 'w', 'Visible', 'off'); 
    set(fig, 'Position', [100, 100, 750, 550]); 
    hold on; grid on;
    set(gca, 'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.5);
    color_main = [0 0.4470 0.7410];       
    color_highlight = [0.8500 0.3250 0.0980]; 

    plot(TrainT_mean, SMSE_mean, '-o', 'Color', color_main, 'LineWidth', 2.2, ...
        'MarkerSize', 8, 'MarkerFaceColor', color_main, 'MarkerEdgeColor', 'w');

    highlight_idx = find(M_list == target_M);
    if ~isempty(highlight_idx)
        plot(TrainT_mean(highlight_idx), SMSE_mean(highlight_idx), 'p', ... 
            'Color', color_highlight, 'LineWidth', 1.5, ...
            'MarkerSize', 16, 'MarkerFaceColor', color_highlight);
    end

    x_range = max(TrainT_mean) - min(TrainT_mean);
    x_offset = x_range * 0.02; 

    for mi = 1:numel(M_list)
        if mi == highlight_idx
            text(TrainT_mean(mi) + x_offset, SMSE_mean(mi), sprintf('M=%d (Chosen)', M_list(mi)), ...
                'FontSize', 12, 'Color', color_highlight, 'FontWeight', 'bold');
        else
            text(TrainT_mean(mi) + x_offset, SMSE_mean(mi), sprintf('M=%d', M_list(mi)), ...
                'FontSize', 11, 'Color', 'k');
        end
    end

    xlabel('Train Time (ms/pt)', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('SMSE', 'FontSize', 13, 'FontWeight', 'bold');
    title(sprintf('%s / %s: Train-Time vs SMSE Trade-off', DatasetName, upper(Method)), 'FontSize', 15, 'FontWeight', 'bold');
    set(gca, 'FontSize', 12, 'LineWidth', 1.2);

    exportgraphics(fig, sprintf('Fig_%s_%s_tradeoff.pdf', DatasetName, upper(Method)), 'ContentType', 'vector');
    exportgraphics(fig, sprintf('Fig_%s_%s_tradeoff.png', DatasetName, upper(Method)), 'Resolution', 300);
    try saveas(fig, sprintf('Fig_%s_%s_tradeoff.emf', DatasetName, upper(Method))); catch; end
    close(fig); 
end

% =========================================================================
% 🖨️ 最终时刻：打印完美的 Markdown 排版代码
% =========================================================================
fprintf('\n\n=======================================================================\n');
fprintf('👇 请复制以下所有内容，粘贴到你的 Markdown 编辑器 (如 Typora) 中 👇\n');
fprintf('=======================================================================\n\n');

fprintf('### Table 1: Ablation Study on Inducing Points (M)\n\n');

% 🌟 终极占位符大法：使用不可见的空格撑大后三列，从而压缩第一列
fprintf('| Dataset | Inducing<br>Points | Train Time (ms/pt)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Test Time (ms/pt)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | SMSE&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; |\n');
fprintf('| :--- | :--- | :--- | :--- | :--- |\n');
for i = 1:length(all_table_rows)
    fprintf('%s\n', all_table_rows{i});
end

fprintf('\n---\n\n');

fprintf('### Figure X: Computational Cost vs. Performance Trade-off\n\n');
fprintf('| | |\n');
fprintf('| :---: | :---: |\n');
fprintf('| <img src="Fig_KIN40K_POE_tradeoff.png" width="450"> | <img src="Fig_POL_POE_tradeoff.png" width="450"> |\n');
fprintf('| <img src="Fig_PUMADYN32NM_POE_tradeoff.png" width="450"> | <img src="Fig_SARCOS_POE_tradeoff.png" width="450"> |\n\n');

fprintf('=======================================================================\n\n');