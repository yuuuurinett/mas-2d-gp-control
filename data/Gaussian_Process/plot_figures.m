% plot_all_curves_final.m


clc; clear; close all;

datasets = {'KIN40K', 'POL', 'PUMADYN32NM', 'SARCOS'};
metrics  = {'SMSE', 'RMSE'};
train_ratio = 40;
n_mc = 5;

% 定义 5 张图的配置
fig_configs = {
    'Fig1_IP_AC',   {'LoG-MOE (Baseline)', 'MOE','POE','GPOE','BCM','RBCM'}, {'log_moe_mc%d.mat', 'moe_ac_tr%d_mc%d.mat','poe_ac_tr%d_mc%d.mat','gpoe_ac_tr%d_mc%d.mat','bcm_ac_tr%d_mc%d.mat','rbcm_ac_tr%d_mc%d.mat'};
    'Fig2_IP_DAC',  {'LoG-MOE (Baseline)', 'MOE','POE','GPOE','BCM','RBCM'}, {'log_moe_mc%d.mat', 'moe_tr%d_mc%d.mat','poe_tr%d_mc%d.mat','gpoe_tr%d_mc%d.mat','bcm_tr%d_mc%d.mat','rbcm_tr%d_mc%d.mat'};
    'Fig3_TP_AC',   {'LoG-MOE (Baseline)', 'MOE','POE','GPOE','BCM','RBCM'}, {'log_moe_mc%d.mat', 'moe_ac_tp_tr%d_mc%d.mat','poe_ac_tp_tr%d_mc%d.mat','gpoe_ac_tp_tr%d_mc%d.mat','bcm_ac_tp_tr%d_mc%d.mat','rbcm_ac_tp_tr%d_mc%d.mat'};
    'Fig4_TP_DAC',  {'LoG-MOE (Baseline)', 'MOE','POE','GPOE','BCM','RBCM'}, {'log_moe_mc%d.mat', 'moe_tp_tr%d_mc%d.mat','poe_tp_tr%d_mc%d.mat','gpoe_tp_tr%d_mc%d.mat','bcm_tp_tr%d_mc%d.mat','rbcm_tp_tr%d_mc%d.mat'};
    'Fig5_Neighbor',{'LoG-MOE (Baseline)', 'MOE','POE','GPOE','BCM','RBCM'}, {'log_moe_mc%d.mat', 'moe_nbr_mc%d.mat','poe_nbr_mc%d.mat','gpoe_nbr_mc%d.mat','bcm_nbr_mc%d.mat','rbcm_nbr_mc%d.mat'};
};

colors = lines(7); 

for d = 1:numel(datasets)
    ds = datasets{d};
    SaveFolder = fullfile('Result', 'Dataset', ds);
    if ~exist(SaveFolder, 'dir'), mkdir(SaveFolder); end
    fprintf('\n========== 正在绘制数据集: %s ==========\n', ds);
    
    for metric_idx = 1:numel(metrics)
        metric = metrics{metric_idx};
        curve_var = lower([metric, '_curve']); 
        mean_var  = lower(metric);             
        
        for f = 1:size(fig_configs, 1)
            fig_name = fig_configs{f, 1};
            labels   = fig_configs{f, 2};
            files    = fig_configs{f, 3};
            
            fig = figure('Name', sprintf('%s_%s_%s', ds, fig_name, metric), 'Position', [100, 100, 900, 600], 'Color', 'w');
            ax_main = axes('Parent', fig, 'Position', [0.10, 0.12, 0.86, 0.80]);
            hold(ax_main, 'on'); grid(ax_main, 'on');
            
            valid_plot = false; max_x_length = 0;   
            
            for m_idx = 1:numel(labels)
                fname_template = files{m_idx};
                all_mc_curves = [];
                for mc = 1:n_mc
                    if contains(fname_template, 'log_') || contains(fname_template, '_nbr_')
                        fname = fullfile(SaveFolder, sprintf(fname_template, mc));
                    else
                        fname = fullfile(SaveFolder, sprintf(fname_template, train_ratio, mc));
                    end
                    if exist(fname, 'file')
                        data = load(fname);
                        if isfield(data, curve_var)
                            cur_c = data.(curve_var);
                            if size(cur_c, 1) > size(cur_c, 2), cur_c = cur_c'; end 
                            all_mc_curves = [all_mc_curves; cur_c];
                        end
                    end
                end
                
                if ~isempty(all_mc_curves)
                    mean_curve = mean(all_mc_curves, 1, 'omitnan');
                    x_axis = 1:length(mean_curve);
                    max_x_length = max(max_x_length, length(mean_curve));
                    line_style = '-'; line_width = 2.0; 
                    if contains(labels{m_idx}, 'Baseline')
                        line_color = [0.6 0.6 0.6]; line_style = '--'; line_width = 3.0; 
                    else
                        line_color = colors(m_idx - 1, :);
                        if contains(labels{m_idx}, 'RBCM'), line_style = '-.'; line_width = 2.5; end
                    end
                    plot(ax_main, x_axis, mean_curve, 'LineStyle', line_style, 'LineWidth', line_width, 'Color', line_color, 'DisplayName', labels{m_idx});
                    valid_plot = true;
                end
            end
            
            if valid_plot
               
                main_lines = findobj(ax_main, 'Type', 'line');
                zoom_x_start = round(max_x_length * 0.7); zoom_x_end = max_x_length;
                focus_y_min = inf; focus_y_max = -inf; all_y_max = -inf;
                
                for i = 1:length(main_lines)
                    ydata = main_lines(i).YData; xdata = main_lines(i).XData; name = main_lines(i).DisplayName;
                    valid_idx = ~isnan(ydata) & ~isinf(ydata);
                    if ~any(valid_idx), continue; end
                    idx_stable = (xdata > max_x_length * 0.05) & valid_idx;
                    if any(idx_stable), all_y_max = max(all_y_max, max(ydata(idx_stable))); end
                    idx_zoom = (xdata >= zoom_x_start) & valid_idx;
                    if any(idx_zoom) && ~contains(name, 'Baseline')
                        focus_y_min = min(focus_y_min, min(ydata(idx_zoom)));
                        focus_y_max = max(focus_y_max, max(ydata(idx_zoom)));
                    end
                end
                

                y_span = all_y_max - focus_y_min;
                if isinf(y_span) || isnan(y_span) || y_span <= 0, y_span = focus_y_min * 0.1; end
                ylim(ax_main, [focus_y_min - y_span * 1.8, all_y_max + y_span * 0.15]);
                y_ticks = get(ax_main, 'YTick');
                set(ax_main, 'YTick', y_ticks(y_ticks >= focus_y_min - y_span*0.05));
                
               
                title(ax_main, sprintf('%s - %s - %s', ds, strrep(fig_name, '_', ' '), metric), 'FontSize', 16, 'FontWeight', 'bold');
                xlabel(ax_main, 'Number of Evaluated Test Points', 'FontSize', 14, 'FontWeight', 'bold');
                ylabel(ax_main, metric, 'FontSize', 14, 'FontWeight', 'bold');
                legend(ax_main, 'Location', 'northeast', 'FontSize', 11, 'NumColumns', 2); 
                ax_main.XAxis.Exponent = 0;
                
               
                if focus_y_max > -inf
                    focus_span = focus_y_max - focus_y_min;
                    y_min_zoom = focus_y_min - focus_span*0.15; y_max_zoom = focus_y_max + focus_span*0.15;
                    patch(ax_main, 'XData', [zoom_x_start, zoom_x_end, zoom_x_end, zoom_x_start], ...
                          'YData', [y_min_zoom, y_min_zoom, y_max_zoom, y_max_zoom], ...
                          'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                    
                    ax_inset = axes('Parent', fig, 'Position', [0.18, 0.25, 0.70, 0.32]); 
                    box(ax_inset, 'on'); hold(ax_inset, 'on'); grid(ax_inset, 'on');
                    copyobj(main_lines, ax_inset);
                    xlim(ax_inset, [zoom_x_start, zoom_x_end]); ylim(ax_inset, [y_min_zoom, y_max_zoom]);
                    set(ax_inset, 'FontSize', 11, 'LineWidth', 1.2, 'Color', 'w');
                    ax_inset.XAxis.Exponent = 0;
                end

               
                base_name = fullfile(SaveFolder, sprintf('%s_%s', fig_name, metric));
                % 1. 导出高清 PNG (推荐用于快速预览和插入 PPT)
                exportgraphics(fig, [base_name, '.png'], 'Resolution', 600);
                % 2. 导出矢量 EMF (Windows/PPT 官方推荐，无限清晰)
                exportgraphics(fig, [base_name, '.emf'], 'ContentType', 'vector');
                
                fprintf('  已生成: %s.png/emf\n', sprintf('%s_%s', fig_name, metric));
            end
            close(fig);
        end
    end
end
