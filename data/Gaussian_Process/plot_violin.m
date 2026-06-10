
% plot_ultimate_right_inset_perfect.m
clc; clear; close all;

% =========================================================================
% 🎛️ 绝对精准的微观刻度
% =========================================================================
% 依次对应: KIN40K, POL, PUMADYN32NM, SARCOS
INSET_MAX_TRAIN = [0.12, 0.35, 0.35, 2.45]; 
INSET_MAX_TEST  = [0.15, 0.40, 0.35, 0.60]; 

datasets = {'KIN40K', 'POL', 'PUMADYN32NM', 'SARCOS'};
aggs     = {'moe', 'gpoe', 'poe', 'bcm', 'rbcm'};
train_ratio = 40;
n_mc = 1; 

archs = {
    'LoG',    'log_%s_mc%d.mat';
    'CEN',    '%s_cen_tr%d_mc%d.mat';  
    'IP-DAC', '%s_tr%d_mc%d.mat';
    'IP-AC',  '%s_ac_tr%d_mc%d.mat';
    'TP-DAC', '%s_tp_tr%d_mc%d.mat';
    'TP-AC',  '%s_ac_tp_tr%d_mc%d.mat';
    'NBR',    '%s_nbr_tr%d_mc%d.mat'
};
method_names = archs(:, 1)';

colors_7 = [
    0.65, 0.65, 0.65; 0.35, 0.35, 0.35; 0.85, 0.33, 0.10;
    0.00, 0.45, 0.74; 0.93, 0.69, 0.13; 0.49, 0.18, 0.56; 0.46, 0.67, 0.18
];

% ================= 数据读取 =================
all_train_data = cell(4, 7); all_test_data  = cell(4, 7);
fprintf('开始读取数据...\n');
for d_idx = 1:numel(datasets)
    ds = datasets{d_idx}; SaveFolder = fullfile('Result', 'Dataset', ds);
    for arch_idx = 1:size(archs, 1)
        fname_tpl = archs{arch_idx, 2}; tr_temp = []; te_temp = [];
        for a_idx = 1:numel(aggs)
            agg = aggs{a_idx};
            for mc = 1:n_mc
                if contains(fname_tpl, 'tr%d'), target_f = fullfile(SaveFolder, sprintf(fname_tpl, agg, train_ratio, mc));
                else, target_f = fullfile(SaveFolder, sprintf(fname_tpl, agg, mc)); end
                
                if exist(target_f, 'file')
                    data = load(target_f); t_tr = []; t_te = [];
                    if isfield(data, 't_train_per_point'), t_tr = data.t_train_per_point;
                    elseif isfield(data, 't_train'),       t_tr = data.t_train; end
                    if isfield(data, 't_test_per_point'),  t_te = data.t_test_per_point;
                    elseif isfield(data, 't_test'),        t_te = data.t_test; end
                    
                    if ~isempty(t_tr) && t_tr > 50 && ~isfield(data, 't_train_per_point')
                        if strcmpi(ds,'KIN40K'), t_tr = (t_tr/40000)*1000; else, t_tr = t_tr*0.01; end
                    end
                    if ~isempty(t_te) && t_te > 50 && ~isfield(data, 't_test_per_point')
                        t_te = (t_te/3000)*1000; 
                    end
                    if ~isempty(t_tr), tr_temp(end+1) = t_tr; end
                    if ~isempty(t_te), te_temp(end+1) = t_te; end
                end
            end
        end
        all_train_data{d_idx, arch_idx} = tr_temp; all_test_data{d_idx, arch_idx}  = te_temp;
    end
end

% ================= 绘图逻辑 =================
for d_idx = 1:numel(datasets)
    ds_name = datasets{d_idx};
    fig = figure('Name', sprintf('%s Results', ds_name), 'Position', [100, 100, 1600, 750], 'Color', 'w');
    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    %% ----------------- 左侧：Train Time -----------------
    ax_t = nexttile(t); hold(ax_t, 'on'); grid(ax_t, 'on');
    all_tr_y = [all_train_data{d_idx, :}]; all_tr_y = all_tr_y(~isnan(all_tr_y));
    inset_max_tr = INSET_MAX_TRAIN(d_idx);
    
    if ~isempty(all_tr_y)
        y_max_t = max(all_tr_y); 
        % 舒适的2.2倍拉升，完美将画中画托起在半空中
        ylim(ax_t, [-0.02, y_max_t * 2.20]); 
    else, y_max_t = 1; end
    
    % 主图提示框（精准只框 CEN 到 NBR）
    plot(ax_t, [1.5, 7.8, 7.8, 1.5, 1.5], [-0.01, -0.01, inset_max_tr, inset_max_tr, -0.01], 'k--', 'LineWidth', 1.2);
    
    for i = 1:7, draw_core_violin(ax_t, i, all_train_data{d_idx, i}, colors_7(i, :), false, y_max_t); end
    
    set(ax_t, 'FontSize', 15, 'XTick', 1:7, 'XTickLabel', method_names, 'LineWidth', 1.5);
    xtickangle(ax_t, 25); ylabel(ax_t, 'Train Time (ms/pt)', 'FontSize', 18, 'FontWeight', 'bold');
    
    title(ax_t, 'Train Time', 'FontSize', 22, 'FontWeight', 'bold');
    xlim(ax_t, [0.2, 7.8]);
    
    % --- Train 画中画 ---
    drawnow; ax_pos_t = ax_t.Position;
    % 【数学级对齐】精准计算宽度，保证画中画的X轴刻度与下方主图严格垂直对齐
    inset_left_t  = ax_pos_t(1) + ax_pos_t(3) * 0.1710;
    inset_width_t = ax_pos_t(3) * 0.8290;
    inset_ax_t = axes('Position', [inset_left_t, ax_pos_t(2)+ax_pos_t(4)*0.55, inset_width_t, ax_pos_t(4)*0.35]);
    hold(inset_ax_t, 'on'); box(inset_ax_t, 'on'); grid(inset_ax_t, 'on'); set(inset_ax_t, 'Color', 'w');
    
    for i = 1:7, draw_core_violin(inset_ax_t, i, all_train_data{d_idx, i}, colors_7(i, :), true, inset_max_tr); end
    
    ylim(inset_ax_t, [-0.005, inset_max_tr]); xlim(inset_ax_t, [1.5, 7.8]);
    set(inset_ax_t, 'FontSize', 13, 'XTick', 2:7, 'XTickLabel', {}, 'LineWidth', 1.2, 'Color', [0.99 0.99 0.99]);

    %% ----------------- 右侧：Test Time -----------------
    ax_e = nexttile(t); hold(ax_e, 'on'); grid(ax_e, 'on');
    all_te_y = [all_test_data{d_idx, :}]; all_te_y = all_te_y(~isnan(all_te_y));
    inset_max_te = INSET_MAX_TEST(d_idx);
    
    if ~isempty(all_te_y)
        y_max_e = max(all_te_y); 
        ylim(ax_e, [-0.02, y_max_e * 2.20]); 
    else, y_max_e = 1; end
    
    plot(ax_e, [1.5, 7.8, 7.8, 1.5, 1.5], [-0.01, -0.01, inset_max_te, inset_max_te, -0.01], 'k--', 'LineWidth', 1.2);
    
    for i = 1:7, draw_core_violin(ax_e, i, all_test_data{d_idx, i}, colors_7(i, :), false, y_max_e); end
    
    set(ax_e, 'FontSize', 15, 'XTick', 1:7, 'XTickLabel', method_names, 'LineWidth', 1.5);
    xtickangle(ax_e, 25); ylabel(ax_e, 'Test Time (ms/pt)', 'FontSize', 18, 'FontWeight', 'bold');
    
    title(ax_e, 'Test Time', 'FontSize', 22, 'FontWeight', 'bold');
    xlim(ax_e, [0.2, 7.8]);
    
    % --- Test 画中画 ---
    drawnow; ax_pos_e = ax_e.Position;
    inset_left_e  = ax_pos_e(1) + ax_pos_e(3) * 0.1710;
    inset_width_e = ax_pos_e(3) * 0.8290;
    inset_ax_e = axes('Position', [inset_left_e, ax_pos_e(2)+ax_pos_e(4)*0.55, inset_width_e, ax_pos_e(4)*0.35]);
    hold(inset_ax_e, 'on'); box(inset_ax_e, 'on'); grid(inset_ax_e, 'on'); set(inset_ax_e, 'Color', 'w');
    
    for i = 1:7, draw_core_violin(inset_ax_e, i, all_test_data{d_idx, i}, colors_7(i, :), true, inset_max_te); end
    
    ylim(inset_ax_e, [-0.005, inset_max_te]); xlim(inset_ax_e, [1.5, 7.8]);
    set(inset_ax_e, 'FontSize', 13, 'XTick', 2:7, 'XTickLabel', {}, 'LineWidth', 1.2, 'Color', [0.99 0.99 0.99]);

    % 唯一总标题
    sgtitle(fig, sprintf('Computational Cost: %s Dataset', ds_name), 'FontSize', 26, 'FontWeight', 'bold');
    
    saveas(fig, sprintf('Fig_%s_Final_Perfect.png', ds_name));
end

% =========================================================================
% 🛡️ 纯净内核渲染逻辑 (100%保留数据，不再做删减)
% =========================================================================
function draw_core_violin(ax, x_pos, Y, clr, is_inset, current_max)
    % 仅在画中画里跳过 LoG
    if is_inset && x_pos == 1
        return; 
    end

    if isempty(Y) || all(isnan(Y))
        scatter(ax, x_pos, 0.01, 30, 'o', 'MarkerEdgeColor', [0.8 0.8 0.8], 'LineWidth', 1.0);
        return; 
    end
    Y = Y(:); med = median(Y);
    
    if is_inset, violin_w = 0.48; dot_sz = 60; sc_sz = 20; lw_box = 4.5;
    else,        violin_w = 0.45; dot_sz = 70; sc_sz = 20; lw_box = 7.0; end

    if max(Y) == min(Y) || length(Y) < 3
        bw = current_max * 0.02;
        if is_inset, bw = current_max * 0.06; end 
        xi = linspace(med - bw*2.5, med + bw*2.5, 100);
        f = normpdf(xi, med, bw*0.8);
        f = f / max(f) * violin_w; 
        patch(ax, [x_pos + f, fliplr(x_pos - f)], [xi, fliplr(xi)], clr, 'EdgeColor', 'none', 'FaceAlpha', 0.65);
        plot(ax, [x_pos, x_pos], [med-bw*1.5, med+bw*1.5], 'k-', 'LineWidth', 1.2); 
        plot(ax, [x_pos, x_pos], [med-bw*0.5, med+bw*0.5], 'k-', 'LineWidth', lw_box); 
        scatter(ax, x_pos, med, dot_sz, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        return;
    end
    
    Y_jittered = Y + randn(size(Y)) * 1e-6; 
    if is_inset, bw = max(std(Y_jittered)/1.2, current_max * 0.06); 
    else,        bw = max(std(Y_jittered)/1.2, current_max * 0.015); end
    
    [f, xi] = ksdensity(Y_jittered, 'Bandwidth', bw);
    valid_mask = (xi >= min(Y) - bw*1.5) & (xi <= max(Y) + bw*1.5) & (xi >= 0);
    f = f(valid_mask); xi = xi(valid_mask);
    if isempty(f), return; end
    f = f / max(f) * violin_w;
    
    patch(ax, [x_pos + f, fliplr(x_pos - f)], [xi, fliplr(xi)], clr, 'EdgeColor', 'none', 'FaceAlpha', 0.65);
    
    q1 = prctile(Y, 25); q3 = prctile(Y, 75);
    plot(ax, [x_pos, x_pos], [min(Y), max(Y)], 'k-', 'LineWidth', 1.5); 
    plot(ax, [x_pos, x_pos], [q1, q3], 'k-', 'LineWidth', lw_box); 
    
    x_scatter = x_pos + (rand(size(Y)) - 0.5) * 0.15; 
    scatter(ax, x_scatter, Y, sc_sz, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
    scatter(ax, x_pos, med, dot_sz, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5); 
end