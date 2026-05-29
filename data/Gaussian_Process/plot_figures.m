% plot_architectures_comparison_curves.m
%clc; clear; close all;

%datasets = {'KIN40K', 'POL', 'PUMADYN32NM', 'SARCOS'};
datasets = {'SARCOS'};
metrics  = {'SMSE', 'RMSE'};
aggs     = {'moe', 'gpoe', 'poe', 'bcm', 'rbcm'};
train_ratio = 40;
n_mc = 5; % 根据你实际跑的 Monte Carlo 次数修改

% =========================================================================
% [终极排版配置字典]
% 格式: {'数据集', [画中画 X, Y, 宽, 高], [放大起始比例, 结束比例], '图例位置'}
% =========================================================================
dataset_configs = {
    % 格式: '数据集', [画中画 X, Y, 宽, 高], [放大起始比例, 结束比例], '图例位置'
    'KIN40K',      [0.20, 0.18, 0.35, 0.35], [0.01, 0.15], 'northeast';

    'POL',         [0.22, 0.22, 0.42, 0.38], [0.02, 0.10], 'northeast';

    'PUMADYN32NM', [0.20, 0.18, 0.38, 0.35], [0.03, 0.18], 'northeast';

    'SARCOS',      [0.55, 0.12, 0.38, 0.32], [0.70, 1.00], 'northwest';
    };

arch_configs = {
    'LoG (Baseline)', 'log_%s_mc%d.mat',        false, '-',  3.0;
    'CEN',            'cen_%s_tr%d_mc%d.mat',   true,  '-',  2.5; 
    'IP-DAC (Ours)',  '%s_tr%d_mc%d.mat',       true,  '-',  2.5; 
    'IP-AC',          '%s_ac_tr%d_mc%d.mat',    true,  '--', 2.0;
    'TP-DAC',         '%s_tp_tr%d_mc%d.mat',    true,  '-.', 2.0;
    'TP-AC',          '%s_ac_tp_tr%d_mc%d.mat', true,  ':',  2.0;
    'NBR',            '%s_nbr_tr%d_mc%d.mat',   true,  '-',  2.0; 
};

colors = [
    0.60, 0.60, 0.60; % LoG: 灰色
    0.00, 0.00, 0.00; % CEN: 黑色
    0.85, 0.33, 0.10; % IP-DAC: 亮橙/红色 (突出显示 Ours)
    0.00, 0.45, 0.74; % IP-AC: 经典蓝
    0.93, 0.69, 0.13; % TP-DAC: 金色
    0.49, 0.18, 0.56; % TP-AC: 紫色
    0.47, 0.67, 0.19  % NBR: 绿色
];

for d = 1:numel(datasets)
    ds = datasets{d};
    SaveFolder = fullfile('Result', 'Dataset', ds);
    if ~exist(SaveFolder, 'dir'), continue; end
    
    ds_idx = find(strcmp(dataset_configs(:,1), ds));
    if isempty(ds_idx), continue; end
    inset_pos  = dataset_configs{ds_idx, 2};
    zoom_range = dataset_configs{ds_idx, 3}; % 现在是一个数组 [start, end]
    legend_loc = dataset_configs{ds_idx, 4};
    
    fprintf('\n========== 正在绘制数据集: %s ==========\n', ds);
    
    for metric_idx = 1:numel(metrics)
        metric = metrics{metric_idx};
        curve_var = lower([metric, '_curve']); 
        
        for a_idx = 1:numel(aggs)
            agg = aggs{a_idx};
            
            fig_name = sprintf('%s_%s_%s', ds, upper(agg), metric);
            fig = figure('Name', fig_name, 'Position', [100, 100, 950, 650], 'Color', 'w');
            ax_main = axes('Parent', fig, 'Position', [0.10, 0.12, 0.86, 0.80]);
            hold(ax_main, 'on'); grid(ax_main, 'on');
            
            curves_buffer = cell(size(arch_configs, 1), n_mc);
            min_curve_len = inf;
            
            for arch_idx = 1:size(arch_configs, 1)
                fname_tpl = arch_configs{arch_idx, 2};
                uses_tr   = arch_configs{arch_idx, 3};
                for mc = 1:n_mc
                    if uses_tr, fname = fullfile(SaveFolder, sprintf(fname_tpl, agg, train_ratio, mc));
                    else,       fname = fullfile(SaveFolder, sprintf(fname_tpl, agg, mc)); end
                    
                    if exist(fname, 'file')
                        data = load(fname);
                        if isfield(data, curve_var)
                            cur_c = data.(curve_var);
                            if size(cur_c, 1) > size(cur_c, 2), cur_c = cur_c'; end 
                            curves_buffer{arch_idx, mc} = cur_c;
                            min_curve_len = min(min_curve_len, length(cur_c));
                        end
                    end
                end
            end
            
            if min_curve_len == inf 
                close(fig); continue; 
            end 
            
            valid_plot = false;
            max_x_length = min_curve_len; 
            
            for arch_idx = 1:size(arch_configs, 1)
                arch_name = arch_configs{arch_idx, 1};
                line_sty  = arch_configs{arch_idx, 4};
                line_wid  = arch_configs{arch_idx, 5};
                
                all_mc_curves = [];
                for mc = 1:n_mc
                    if ~isempty(curves_buffer{arch_idx, mc})
                        cur_c = curves_buffer{arch_idx, mc};
                        all_mc_curves = [all_mc_curves; cur_c(1:min_curve_len)]; 
                    end
                end
                
                if ~isempty(all_mc_curves)
                    mean_curve = mean(all_mc_curves, 1, 'omitnan');
                    x_axis = 1:max_x_length;
                    
                    plot(ax_main, x_axis, mean_curve, 'LineStyle', line_sty, ...
                        'LineWidth', line_wid, 'Color', colors(arch_idx, :), ...
                        'DisplayName', arch_name);
                    valid_plot = true;
                end
            end
            
            if valid_plot
                main_lines = findobj(ax_main, 'Type', 'line');
                
                % [核心改动] 精准截取自定义的放大区域
                zoom_x_start = max(1, round(max_x_length * zoom_range(1))); 
                zoom_x_end   = min(max_x_length, round(max_x_length * zoom_range(2)));
                
                focus_y_min = inf; focus_y_max = -inf; all_y_max = -inf;
                
                for i = 1:length(main_lines)
                    ydata = main_lines(i).YData; xdata = main_lines(i).XData; name = main_lines(i).DisplayName;
                    valid_idx = ~isnan(ydata) & ~isinf(ydata);
                    if ~any(valid_idx), continue; end
                    
                    % 主图 Y 轴上限计算
                    idx_stable = (xdata > max_x_length * 0.05) & valid_idx;
                    if any(idx_stable), all_y_max = max(all_y_max, max(ydata(idx_stable))); end
                    
                    % -------------------------------------------------
                    % [高级优化] 仅在自定义的 zoom 区间内寻找放大框的 Y 轴极值
                    % 核心技巧：计算放大框 Y 轴时，主动忽略高高在上的 NBR 和 Baseline(LoG)
                    % 这样放大框的 Y 轴就会自动缩窄，把中间和底部的 IP/TP 曲线纵向放大拉开！
                    % -------------------------------------------------
                    idx_zoom = (xdata >= zoom_x_start & xdata <= zoom_x_end) & valid_idx;
                    if any(idx_zoom) && ~contains(name, 'LoG') && ~contains(name, 'NBR')
                        focus_y_min = min(focus_y_min, min(ydata(idx_zoom)));
                        focus_y_max = max(focus_y_max, max(ydata(idx_zoom)));
                    end
                end
                y_span = all_y_max - focus_y_min;
                if isinf(y_span) || isnan(y_span) || y_span <= 0, y_span = focus_y_min * 0.1; end
                ylim(ax_main, [focus_y_min - y_span * 1.5, all_y_max + y_span * 0.15]);
                
                title(ax_main, sprintf('%s Dataset - %s Aggregation - %s Curve', ds, upper(agg), metric), ...
                    'FontSize', 16, 'FontWeight', 'bold');
                xlabel(ax_main, 'Number of Evaluated Test Points', 'FontSize', 14, 'FontWeight', 'bold');
                ylabel(ax_main, metric, 'FontSize', 14, 'FontWeight', 'bold');
                legend(ax_main, 'Location', legend_loc, 'FontSize', 11, 'NumColumns', 2); 
                ax_main.XAxis.Exponent = 0;
                
                if focus_y_max > -inf && focus_y_min < inf
                    focus_span = focus_y_max - focus_y_min;
                    if focus_span == 0, focus_span = focus_y_max * 0.01; end
                    y_min_zoom = focus_y_min - focus_span * 0.2; 
                    y_max_zoom = focus_y_max + focus_span * 0.2;
                    
                    patch(ax_main, 'XData', [zoom_x_start, zoom_x_end, zoom_x_end, zoom_x_start], ...
                          'YData', [y_min_zoom, y_min_zoom, y_max_zoom, y_max_zoom], ...
                          'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                    
                    ax_inset = axes('Parent', fig, 'Position', inset_pos); 
                    box(ax_inset, 'on'); hold(ax_inset, 'on'); grid(ax_inset, 'on');
                    copyobj(main_lines, ax_inset);
                    xlim(ax_inset, [zoom_x_start, zoom_x_end]); 
                    ylim(ax_inset, [y_min_zoom, y_max_zoom]);
                    set(ax_inset, 'FontSize', 11, 'LineWidth', 1.2, 'Color', 'w');
                    ax_inset.XAxis.Exponent = 0;
                end

                OutName = fullfile(SaveFolder, sprintf('Curve_%s_%s_%s', ds, upper(agg), metric));
                savefig(fig, [OutName, '.fig']);
                fprintf('  已生成并保留窗口: %s.fig\n', sprintf('Curve_%s_%s_%s', ds, upper(agg), metric));
            end
        end
    end
end