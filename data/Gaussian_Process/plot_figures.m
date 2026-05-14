%% ================= 终极批处理：全自动生成 6 张高级阴影对比图 =================
clc; clear; close all;

% 全局参数配置
dataset_name = 'KIN40K';
folder = fullfile('Result', 'Dataset', dataset_name);
N_MC = 5;              % 你实际跑了多少次MC
downsample_rate = 20;  % 降采样率(避免矢量图卡死)

% 配色库 (保持统一，比如 MOE永远是红色，gPOE永远是蓝色)
% 1橙红(MOE), 2深蓝(gPOE), 3黄(POE), 4紫(BCM), 5绿(rBCM), 6深灰(LoG)
colors = [
    0.85, 0.33, 0.10; 
    0.00, 0.45, 0.74; 
    0.93, 0.69, 0.13; 
    0.49, 0.18, 0.56; 
    0.47, 0.67, 0.19; 
    0.30, 0.30, 0.30; 
];

%% 1. 定义 6 个绘图任务 (图名 + 包含的曲线)
% 格式: {'图例名称', '文件命名规则', 颜色序号, '线型'}
plot_tasks = {
    % 任务 1: IP-AC 内部对比
    {'Fig1_IP_AC', {
        {'IP-AC (MOE)',  'moe_ac_tr40_mc%d.mat',  1, '-'};
        {'IP-AC (gPOE)', 'gpoe_ac_tr40_mc%d.mat', 2, '-'};
        {'IP-AC (POE)',  'poe_ac_tr40_mc%d.mat',  3, '-'};
        {'IP-AC (BCM)',  'bcm_ac_tr40_mc%d.mat',  4, '-'};
        {'IP-AC (rBCM)', 'rbcm_ac_tr40_mc%d.mat', 5, '-'};
    }};

    % 任务 2: IP-DAC 内部对比
    {'Fig2_IP_DAC', {
        {'IP-DAC (MOE)',  'moe_tr40_mc%d.mat',  1, '--'};
        {'IP-DAC (gPOE)', 'gpoe_tr40_mc%d.mat', 2, '--'};
        {'IP-DAC (POE)',  'poe_tr40_mc%d.mat',  3, '--'};
        {'IP-DAC (BCM)',  'bcm_tr40_mc%d.mat',  4, '--'};
        {'IP-DAC (rBCM)', 'rbcm_tr40_mc%d.mat', 5, '--'};
    }};

    % 任务 3: TP-AC 内部对比
    {'Fig3_TP_AC', {
        {'TP-AC (MOE)',  'moe_ac_tp_tr40_mc%d.mat',  1, '-'};
        {'TP-AC (gPOE)', 'gpoe_ac_tp_tr40_mc%d.mat', 2, '-'};
        {'TP-AC (POE)',  'poe_ac_tp_tr40_mc%d.mat',  3, '-'};
        {'TP-AC (BCM)',  'bcm_ac_tp_tr40_mc%d.mat',  4, '-'};
        {'TP-AC (rBCM)', 'rbcm_ac_tp_tr40_mc%d.mat', 5, '-'};
    }};

    % 任务 4: TP-DAC 内部对比
    {'Fig4_TP_DAC', {
        {'TP-DAC (MOE)',  'moe_tp_tr40_mc%d.mat',  1, '--'};
        {'TP-DAC (gPOE)', 'gpoe_tp_tr40_mc%d.mat', 2, '--'};
        {'TP-DAC (POE)',  'poe_tp_tr40_mc%d.mat',  3, '--'};
        {'TP-DAC (BCM)',  'bcm_tp_tr40_mc%d.mat',  4, '--'};
        {'TP-DAC (rBCM)', 'rbcm_tp_tr40_mc%d.mat', 5, '--'};
    }};

    % 任务 5: Neighbor 静态对比
    {'Fig5_Neighbor', {
        {'NBR (MOE)',  'moe_nbr_mc%d.mat',  1, '-.'};
        {'NBR (gPOE)', 'gpoe_nbr_mc%d.mat', 2, '-.'};
        {'NBR (POE)',  'poe_nbr_mc%d.mat',  3, '-.'};
        {'NBR (BCM)',  'bcm_nbr_mc%d.mat',  4, '-.'};
        {'NBR (rBCM)', 'rbcm_nbr_mc%d.mat', 5, '-.'};
    }};

    % 任务 6: 巅峰对决 (LoG基线 vs 最好的AC vs 最好的DAC)
    % 假设 MOE 是你方法里表现最好的，把它单拎出来和 LoG 打
    {'Fig6_Highlight_Showdown', {
        {'LoG-MOE',      'log_moe_mc%d.mat',      6, '-.'};
        {'LoG-gPOE',     'log_gpoe_mc%d.mat',     6, ':'};
        {'IP-AC (MOE)',  'moe_ac_tr40_mc%d.mat',  1, '-'}; 
        {'IP-DAC (MOE)', 'moe_tr40_mc%d.mat',     1, '--'};
    }};
};

%% 2. 循环执行所有绘图任务
for task_idx = 1:length(plot_tasks)
    fig_name = plot_tasks{task_idx}{1};
    methods_to_plot = plot_tasks{task_idx}{2};
    
    % 初始化画布
    fig = figure('Name', fig_name, 'Position', [100, 100, 650, 450], 'Color', 'w');
    hold on; grid on; box on;
    set(gca, 'YScale', 'log', 'FontSize', 12, 'GridAlpha', 0.25); 
    
    legend_handles = [];
    legend_labels = {};
    
    % 画每一条线
    for i = 1:size(methods_to_plot, 1)
        plot_name = methods_to_plot{i}{1};
        file_tpl  = methods_to_plot{i}{2};
        c_idx     = methods_to_plot{i}{3};
        line_sty  = methods_to_plot{i}{4};
        c_rgb     = colors(c_idx, :);
        
        % 读取数据
        mc_curves = [];
        for mc = 1:N_MC
            fname = fullfile(folder, sprintf(file_tpl, mc));
            if exist(fname, 'file')
                data = load(fname, 'smse_curve');
                if isfield(data, 'smse_curve')
                    mc_curves = [mc_curves, data.smse_curve(:)];
                end
            end
        end
        
        if isempty(mc_curves)
            fprintf('  [跳过] 找不到数据: %s\n', sprintf(file_tpl, 1));
            continue;
        end
        
        % 算均值和方差，并降采样
        mean_curve = mean(mc_curves, 2);
        std_curve  = std(mc_curves, 0, 2);
        
        x_vals = (1 : downsample_rate : length(mean_curve))';
        if x_vals(end) ~= length(mean_curve), x_vals = [x_vals; length(mean_curve)]; end
        
        y_mean = mean_curve(x_vals);
        y_std  = std_curve(x_vals);
        
        % 画半透明阴影
        x_fill = [x_vals; flipud(x_vals)];
        y_fill = [y_mean + y_std; flipud(max(y_mean - y_std, 1e-8))]; 
        fill(x_fill, y_fill, c_rgb, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
        % 画均值实线
        h = plot(x_vals, y_mean, 'LineStyle', line_sty, 'Color', c_rgb, 'LineWidth', 2.0);
        
        legend_handles(end+1) = h;
        legend_labels{end+1}  = plot_name;
    end
    
    % 图表美化
    xlabel('Iteration (Test Points)', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Cumulative SMSE', 'FontSize', 13, 'FontWeight', 'bold');
    title(strrep(fig_name, '_', ' '), 'FontSize', 14, 'FontWeight', 'bold');
    
    xticks = get(gca, 'XTick');
    if max(xticks) >= 1000
        xtick_labels = arrayfun(@(x) sprintf('%dk', x/1000), xticks, 'UniformOutput', false);
        xtick_labels(xticks == 0) = {'0'};
        set(gca, 'XTickLabel', xtick_labels);
    end
    
    if ~isempty(legend_handles)
        lgd = legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 11);
        lgd.Color = [1 1 1 0.8]; 
    end
    
    hold off;
    
    % 输出PDF
    savepath = fullfile(folder, sprintf('%s.pdf', fig_name));
    exportgraphics(fig, savepath, 'ContentType', 'vector');
    fprintf('成功生成: %s\n', savepath);
end
fprintf('\n6 张精美图表已全部生成完毕！快去文件夹里看看吧。\n');