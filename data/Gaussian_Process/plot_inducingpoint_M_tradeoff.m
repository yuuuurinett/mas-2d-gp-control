%% plot_inducingpoint_M_tradeoff_summary.m
% M-ablation 结果汇总图 (参照 Lederer et al. 2021, ICML, Figure 2 风格)
%
% 布局：
%   每个 dataset 一幅图 (共4幅图)
%   每幅图内 上下两排 subplot:
%     上排 = SMSE, log scale, x轴 = M (诱导点数量)
%     下排 = MSLL, 线性 scale, x轴 = M
%   同一排subplot内, 5个method用不同颜色+marker画在一起(legend区分),
%   方便横向对比不同aggregation方法在同一数据集上随M变化的表现
%   数据: mean ± std, 跨 seeds 取统计量, errorbar 展示
clc; close all;

Method_list  = {'poe','gpoe','moe','bcm','rbcm'};
Method_label = {'POE','GPOE','MOE','BCM','RBCM'};
Dataset_list = {'KIN40K','POL','PUMADYN32NM','SARCOS'};
M_list       = 100:100:2500;
seeds        = 1:3;
train_ratio  = 0.4;
tr_tag       = round(train_ratio*100);

ProjectRoot = fileparts(mfilename('fullpath'));

% ---- 1. 读取全部数据，存进结构体 ----
% data.(method).(dataset).SMSE_mean/std, MSLL_mean/std  -- 长度 = numel(M_list)
data = struct();

for ci = 1:numel(Method_list)
    method = Method_list{ci};
    for di = 1:numel(Dataset_list)
        dataset = Dataset_list{di};
        ResultFolder = fullfile(ProjectRoot, 'Result','Dataset',dataset);

        SMSE_all = nan(numel(M_list), numel(seeds));
        MSLL_all = nan(numel(M_list), numel(seeds));

        for mi = 1:numel(M_list)
            M = M_list(mi);
            for si = 1:numel(seeds)
                seed = seeds(si);
                fname = sprintf('%s_M%d_tr%d_mc%d.mat', method, M, tr_tag, seed);
                fpath = fullfile(ResultFolder, fname);
                if ~exist(fpath, 'file')
                    fprintf('[缺失] %s\n', fpath);
                    continue;
                end
                S = load(fpath, 'smse', 'msll');
                if isfield(S,'smse'), SMSE_all(mi,si) = S.smse; end
                if isfield(S,'msll'), MSLL_all(mi,si) = S.msll; end
            end
        end

        data.(method).(dataset).SMSE_mean = mean(SMSE_all, 2, 'omitnan');
        data.(method).(dataset).SMSE_std  = std(SMSE_all, 0, 2, 'omitnan');
        data.(method).(dataset).MSLL_mean = mean(MSLL_all, 2, 'omitnan');
        data.(method).(dataset).MSLL_std  = std(MSLL_all, 0, 2, 'omitnan');
    end
end

% ---- 2. (SMSE用log scale自动适配范围, MSLL同图对比5个method, 用MATLAB自动y轴缩放) ----

% ---- 3. 画图: 每个 dataset 一幅图, 上排SMSE(log) + 下排MSLL(linear), 5个method同图对比 ----
FigFolder = fullfile(ProjectRoot, 'Result','Figures','M_ablation');
if ~exist(FigFolder, 'dir'), mkdir(FigFolder); end

% 数据本身重叠严重, 颜色区分无意义(被遮挡), 改为全黑+marker形状区分,
% 接近Lederer Fig.2用marker(不靠颜色)辨别曲线的思路
method_colors = repmat([0,0,0], numel(Method_list), 1);  % 全部黑色
method_markers = {'o','^','s','d','v'};
marker_every = 2;  % 每隔几个数据点放一个marker, 避免25个点全标记太密

% 自动识别放大区域: 取M_list后1/3区间(差异通常最小最难区分), 用于inset画中画
inset_idx_start = ceil(2/3 * numel(M_list)) + 1;
inset_M_range = [M_list(inset_idx_start), M_list(end)];

% 先不加偏移, 直接对比5条原始曲线的可分辨程度
offset_pct = 0.00;  % <<< 可调参数: 0=不偏移(先看原始效果); 如仍重叠可调到 0.05~0.10
offset_factor = 1 + offset_pct * (0:numel(Method_list)-1);

for di = 1:numel(Dataset_list)
    dataset = Dataset_list{di};

    fig = figure('Color','w','Position',[60,60,950,800], ...
        'Name', sprintf('%s_SMSE_MSLL_vs_M', dataset));

    % ---- 上排: SMSE, log scale ----
    ax1 = subplot(2,1,1);
    hold(ax1, 'on');
    for ci = 1:numel(Method_list)
        method = Method_list{ci};
        d = data.(method).(dataset);
        % SMSE是正数: 乘法偏移保持相对比例一致, 不扭曲log scale下的形状
        smse_offset = d.SMSE_mean * offset_factor(ci);
        sparse_idx = 1:marker_every:numel(M_list);
        % 完整误差棒+细线 (无marker, 避免25个点全标太密)
        errorbar(ax1, M_list, smse_offset, d.SMSE_std, ...
            '-', 'Color', method_colors(ci,:), 'LineWidth', 1.3, ...
            'CapSize', 2, 'HandleVisibility', 'off');
        % 在稀疏点上叠加marker (用于辨识方法, 不重复画线/误差棒)
        plot(ax1, M_list(sparse_idx), smse_offset(sparse_idx), ...
            'LineStyle', 'none', 'Marker', method_markers{ci}, ...
            'Color', method_colors(ci,:), 'MarkerSize', 6, ...
            'DisplayName', Method_label{ci});
    end
    set(ax1, 'YScale', 'log');
    set(ax1, 'XScale', 'log');
    xticks(ax1, [100, 200, 500, 1000, 2000, 2500]);
    grid(ax1, 'on');
    xlim(ax1, [M_list(1), M_list(end)]);
    % 紧贴真实数据范围设置ylim, 避免log scale默认padding过宽导致曲线"扁"
    all_smse_vals = [];
    for ci = 1:numel(Method_list)
        d_tmp = data.(Method_list{ci}).(dataset);
        all_smse_vals = [all_smse_vals; d_tmp.SMSE_mean * offset_factor(ci) - d_tmp.SMSE_std; ...
                                          d_tmp.SMSE_mean * offset_factor(ci) + d_tmp.SMSE_std]; %#ok<AGROW>
    end
    all_smse_vals = all_smse_vals(all_smse_vals > 0 & isfinite(all_smse_vals));
    if ~isempty(all_smse_vals)
        ylim(ax1, [min(all_smse_vals)*0.9, max(all_smse_vals)*1.1]);
    end
    ylabel(ax1, 'SMSE');
    title(ax1, dataset, 'FontWeight','normal');
    legend(ax1, 'Location','best', 'FontSize',8, 'NumColumns', 5);
    hold(ax1, 'off');

    % ---- SMSE 画中画: 放大M较大区间(差异最小最难区分的部分) ----
    ax1_pos = get(ax1, 'Position');
    ax1_inset = axes('Position', [ax1_pos(1)+ax1_pos(3)*0.45, ax1_pos(2)+ax1_pos(4)*0.45, ...
                                    ax1_pos(3)*0.45, ax1_pos(4)*0.45]);
    hold(ax1_inset, 'on');
    inset_mask = M_list >= inset_M_range(1) & M_list <= inset_M_range(2);
    all_smse_inset = [];
    for ci = 1:numel(Method_list)
        method = Method_list{ci};
        d = data.(method).(dataset);
        smse_offset = d.SMSE_mean * offset_factor(ci);
        plot(ax1_inset, M_list(inset_mask), smse_offset(inset_mask), ...
            '-', 'Marker', method_markers{ci}, 'Color', method_colors(ci,:), ...
            'LineWidth', 1, 'MarkerSize', 4, 'MarkerIndices', 1:2:sum(inset_mask));
        all_smse_inset = [all_smse_inset; smse_offset(inset_mask)]; %#ok<AGROW>
    end
    box(ax1_inset, 'on');
    grid(ax1_inset, 'on');
    set(ax1_inset, 'FontSize', 7);
    xlim(ax1_inset, inset_M_range);
    if ~isempty(all_smse_inset) && range(all_smse_inset) > 0
        pad_inset = 0.05 * range(all_smse_inset);
        ylim(ax1_inset, [min(all_smse_inset)-pad_inset, max(all_smse_inset)+pad_inset]);
    end
    hold(ax1_inset, 'off');

    % ---- 下排: MSLL, linear scale ----
    ax2 = subplot(2,1,2);
    hold(ax2, 'on');
    % MSLL常为负数, 用加法偏移(每个method在原值上加一个小常数), 避免乘法
    % 偏移在负数区间反转方向; 偏移量按该dataset所有method的MSLL范围的1%递增
    all_msll_for_offset = [];
    for ci = 1:numel(Method_list)
        all_msll_for_offset = [all_msll_for_offset; data.(Method_list{ci}).(dataset).MSLL_mean]; %#ok<AGROW>
    end
    msll_range = max(all_msll_for_offset) - min(all_msll_for_offset);
    msll_offset_step = offset_pct * msll_range;
    for ci = 1:numel(Method_list)
        method = Method_list{ci};
        d = data.(method).(dataset);
        msll_offset = d.MSLL_mean + msll_offset_step * (ci - 1);
        sparse_idx = 1:marker_every:numel(M_list);
        errorbar(ax2, M_list, msll_offset, d.MSLL_std, ...
            '-', 'Color', method_colors(ci,:), 'LineWidth', 1.3, ...
            'CapSize', 2, 'HandleVisibility', 'off');
        plot(ax2, M_list(sparse_idx), msll_offset(sparse_idx), ...
            'LineStyle', 'none', 'Marker', method_markers{ci}, ...
            'Color', method_colors(ci,:), 'MarkerSize', 6, ...
            'DisplayName', Method_label{ci});
    end
    set(ax2, 'XScale', 'log');
    xticks(ax2, [100, 200, 500, 1000, 2000, 2500]);
    grid(ax2, 'on');
    xlim(ax2, [M_list(1), M_list(end)]);
    % 紧贴真实数据范围设置ylim
    all_msll_vals = [];
    for ci = 1:numel(Method_list)
        d_tmp = data.(Method_list{ci}).(dataset);
        msll_off_tmp = d_tmp.MSLL_mean + msll_offset_step * (ci - 1);
        all_msll_vals = [all_msll_vals; msll_off_tmp - d_tmp.MSLL_std; msll_off_tmp + d_tmp.MSLL_std]; %#ok<AGROW>
    end
    all_msll_vals = all_msll_vals(isfinite(all_msll_vals));
    if ~isempty(all_msll_vals)
        pad = 0.05 * (max(all_msll_vals) - min(all_msll_vals) + eps);
        ylim(ax2, [min(all_msll_vals)-pad, max(all_msll_vals)+pad]);
    end
    xlabel(ax2, 'Number of inducing points M');
    ylabel(ax2, 'MSLL');
    hold(ax2, 'off');

    % ---- MSLL 画中画: 放大M较大区间 ----
    ax2_pos = get(ax2, 'Position');
    ax2_inset = axes('Position', [ax2_pos(1)+ax2_pos(3)*0.45, ax2_pos(2)+ax2_pos(4)*0.10, ...
                                    ax2_pos(3)*0.45, ax2_pos(4)*0.45]);
    hold(ax2_inset, 'on');
    all_msll_inset = [];
    for ci = 1:numel(Method_list)
        method = Method_list{ci};
        d = data.(method).(dataset);
        msll_offset = d.MSLL_mean + msll_offset_step * (ci - 1);
        plot(ax2_inset, M_list(inset_mask), msll_offset(inset_mask), ...
            '-', 'Marker', method_markers{ci}, 'Color', method_colors(ci,:), ...
            'LineWidth', 1, 'MarkerSize', 4, 'MarkerIndices', 1:2:sum(inset_mask));
        all_msll_inset = [all_msll_inset; msll_offset(inset_mask)]; %#ok<AGROW>
    end
    box(ax2_inset, 'on');
    grid(ax2_inset, 'on');
    set(ax2_inset, 'FontSize', 7);
    xlim(ax2_inset, inset_M_range);
    if ~isempty(all_msll_inset) && range(all_msll_inset) > 0
        pad_inset2 = 0.05 * range(all_msll_inset);
        ylim(ax2_inset, [min(all_msll_inset)-pad_inset2, max(all_msll_inset)+pad_inset2]);
    end
    hold(ax2_inset, 'off');

    if offset_pct > 0
        sgtitle(sprintf('%s: SMSE & MSLL vs M (Train=%.0f%%, mean \\pm std over %d seeds, curves offset %.0f%% for visibility)', ...
            dataset, train_ratio*100, numel(seeds), offset_pct*100));
    else
        sgtitle(sprintf('%s: SMSE & MSLL vs M (Train=%.0f%%, mean \\pm std over %d seeds)', ...
            dataset, train_ratio*100, numel(seeds)));
    end

    saveas(fig, fullfile(FigFolder, sprintf('%s_SMSE_MSLL_vs_M.png', dataset)));
    savefig(fig, fullfile(FigFolder, sprintf('%s_SMSE_MSLL_vs_M.fig', dataset)));
end

fprintf('\n完成，共生成 %d 张图（%d datasets，每张上排SMSE+下排MSLL，5个method同图对比），已保存至 %s\n', ...
    numel(Dataset_list), numel(Dataset_list), FigFolder);

%% ---- 4. 计算"拐点区间": SMSE和MSLL都接近各自最优值时的M交集 ----
% 容忍度: 5% (可调整 tol 改变区间宽度)
tol = 0.05;

fprintf('\n========================================================================\n');
fprintf('  Sweet-spot M region per method/dataset\n');
fprintf('  (M range where both SMSE and MSLL are within %.0f%% of their own optimum)\n', tol*100);
fprintf('========================================================================\n');
fprintf('%-6s %-12s %22s %22s %18s\n', 'Agg','Dataset','SMSE-OK range','MSLL-OK range','Sweet spot (M)');
fprintf('%s\n', repmat('-',1,90));

sweet_spot_results = cell(numel(Method_list)*numel(Dataset_list), 6);
row_idx = 0;

for ci = 1:numel(Method_list)
    method = Method_list{ci};
    for di = 1:numel(Dataset_list)
        dataset = Dataset_list{di};
        d = data.(method).(dataset);

        % SMSE 是正数, 越小越好: 可接受范围 SMSE <= minSMSE*(1+tol)
        minSMSE = min(d.SMSE_mean);
        smse_ok_mask = d.SMSE_mean <= minSMSE * (1 + tol);

        % MSLL 通常是负数, 越小(越负)越好: 可接受范围 MSLL <= minMSLL + tol*|minMSLL|
        % (即不比最优值差超过 |minMSLL| 的 tol 比例; 若 minMSLL 接近 0 则退化处理)
        minMSLL = min(d.MSLL_mean);
        msll_ok_mask = d.MSLL_mean <= minMSLL + tol * abs(minMSLL);

        smse_ok_M = M_list(smse_ok_mask);
        msll_ok_M = M_list(msll_ok_mask);

        smse_range_str = sprintf('[%d, %d]', min(smse_ok_M), max(smse_ok_M));
        msll_range_str = sprintf('[%d, %d]', min(msll_ok_M), max(msll_ok_M));

        % 交集
        sweet_M = intersect(smse_ok_M, msll_ok_M);
        if isempty(sweet_M)
            sweet_str = 'EMPTY (no overlap)';
            sweet_lo = NaN; sweet_hi = NaN;
        else
            sweet_lo = min(sweet_M); sweet_hi = max(sweet_M);
            if sweet_lo == sweet_hi
                sweet_str = sprintf('%d', sweet_lo);
            else
                sweet_str = sprintf('[%d, %d]', sweet_lo, sweet_hi);
            end
        end

        fprintf('%-6s %-12s %22s %22s %18s\n', ...
            Method_label{ci}, dataset, smse_range_str, msll_range_str, sweet_str);

        row_idx = row_idx + 1;
        sweet_spot_results(row_idx,:) = {Method_label{ci}, dataset, smse_range_str, msll_range_str, sweet_str, [sweet_lo, sweet_hi]};
    end
end
fprintf('%s\n', repmat('-',1,90));
fprintf('提示: "EMPTY" 表示该方法/数据集下 SMSE 和 MSLL 在 %.0f%% 容忍度内没有共同的M区间,\n', tol*100);
fprintf('      说明两个指标在M选择上存在明显冲突, 需要人工权衡或放宽容忍度 tol 重新计算。\n');
fprintf('若要在 sweet spot 内进一步精细扫描, 可在该区间用 step=10~20 重新跑\n');
fprintf('      run_inducingpoint_M_ablation.m (改 M_list 即可)。\n');

%% ---- 5. 导出 markdown 表格，方便直接粘贴进 Word ----
md_path = fullfile(FigFolder, 'best_M_summary.md');
fid = fopen(md_path, 'w');
fprintf(fid, '## Sweet-spot M region per method/dataset\n\n');
fprintf(fid, 'Train = %.0f%%, mean over %d seeds (seed = %d~%d). ', ...
    train_ratio*100, numel(seeds), seeds(1), seeds(end));
fprintf(fid, 'M range where both SMSE and MSLL are within %.0f%% of their own optimum (M_list = %d:%d:%d). ', ...
    tol*100, M_list(1), M_list(2)-M_list(1), M_list(end));
fprintf(fid, 'Sweet spot = intersection of the SMSE-OK and MSLL-OK ranges.\n\n');
fprintf(fid, '| Agg | Dataset | SMSE-OK range | MSLL-OK range | Sweet spot (M) |\n');
fprintf(fid, '|---|---|---|---|---|\n');
for r = 1:row_idx
    fprintf(fid, '| %s | %s | %s | %s | %s |\n', ...
        sweet_spot_results{r,1}, sweet_spot_results{r,2}, sweet_spot_results{r,3}, ...
        sweet_spot_results{r,4}, sweet_spot_results{r,5});
end
fprintf(fid, '\n*Note: "EMPTY (no overlap)" means SMSE and MSLL do not share a common acceptable M range at this tolerance; consider relaxing `tol` or treating SMSE/MSLL trade-off explicitly for that case.*\n');
fclose(fid);
fprintf('\nMarkdown 表格已导出到: %s\n', md_path);
fprintf('可直接打开该文件复制内容粘贴进 Word（支持 markdown 表格自动转换的版本）。\n');