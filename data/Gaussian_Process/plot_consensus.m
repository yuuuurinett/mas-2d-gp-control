%% plot_consensus_convergence_summary.m
% 汇总图：4个数据集 x 5个方法的收敛曲线对比
% 每个数据集一个子图，每个方法用不同颜色展示 IP-DAC（实线）和 IP-AC（虚线）
% 曲线图：仅用 seed = plot_seed 作为代表绘制（conv_curve 长度随seed变化，不能逐点平均）
% 汇总表：iters_converge / comm_train / final_disagreement 跨 seed_list 取均值±标准差
clc; close all;

datasets    = {'KIN40K','POL','PUMADYN32NM','SARCOS'};
agg_list    = {'moe','gpoe','poe','bcm','rbcm'};
agg_labels  = {'MOE','GPOE','POE','BCM','RBCM'};
train_ratio = 0.4;
seed_list   = 1:10;   % 跑了10遍 MC 的 seed 列表
plot_seed   = 1;      % 曲线图代表性 seed（保持原脚本行为）
tr_tag      = round(train_ratio*100);

colors = lines(5);  % 5种颜色对应5个method

figure('Color','w','Position',[60,60,1400,950]);

summary_table = {};  % 汇总数值表：均值±标准差（跨 seed_list）

for d_idx = 1:length(datasets)
    dname = datasets{d_idx};
    SaveFolder = fullfile('Result','Dataset',dname);

    subplot(2,2,d_idx);
    hold on;
    legend_entries = {};

    for a_idx = 1:length(agg_list)
        agg = agg_list{a_idx};

        %% ---- 曲线图：只用 plot_seed 作代表，原脚本逻辑不变 ----
        dac_file_plot = fullfile(SaveFolder, sprintf('%s_tr%d_mc%d.mat',    agg, tr_tag, plot_seed));
        ac_file_plot  = fullfile(SaveFolder, sprintf('%s_ac_tr%d_mc%d.mat', agg, tr_tag, plot_seed));
        if exist(dac_file_plot,'file') && exist(ac_file_plot,'file')
            d_dac_plot = load(dac_file_plot);
            d_ac_plot  = load(ac_file_plot);

            iters_dac_plot = d_dac_plot.iter_converge;
            iters_ac_plot  = d_ac_plot.iter_converge;

            % 多条曲线数值上几乎重叠（不同 aggregation 方法在当前实现中
            % consensus 动力学非常接近），加一个很小的相对偏移让线条可分辨，
            % 不影响数量级/收敛趋势的判读。偏移量随 a_idx 递增。
            offset_factor = 1 + 0.05 * (a_idx - 1);

            curve_dac_offset = d_dac_plot.conv_curve_dac * offset_factor;
            curve_ac_offset  = d_ac_plot.conv_curve_ac  * offset_factor;

            semilogy(1:iters_dac_plot, curve_dac_offset, '-',  'Color', colors(a_idx,:), 'LineWidth',1.6);
            semilogy(1:iters_ac_plot,  curve_ac_offset,   '--', 'Color', colors(a_idx,:), 'LineWidth',1.6);

            legend_entries{end+1} = sprintf('%s-DAC', agg_labels{a_idx});
            legend_entries{end+1} = sprintf('%s-AC',  agg_labels{a_idx});
        else
            fprintf('[跳过曲线] %s - %s (seed=%d): 文件不存在\n', dname, agg, plot_seed);
        end

        %% ---- 汇总表：跨 seed_list 收集标量指标，取均值±标准差 ----
        iters_dac_all = []; comm_dac_all = []; final_dac_all = [];
        iters_ac_all  = []; comm_ac_all  = []; final_ac_all  = [];

        for s = seed_list
            dac_file = fullfile(SaveFolder, sprintf('%s_tr%d_mc%d.mat',    agg, tr_tag, s));
            ac_file  = fullfile(SaveFolder, sprintf('%s_ac_tr%d_mc%d.mat', agg, tr_tag, s));

            if ~exist(dac_file,'file') || ~exist(ac_file,'file')
                fprintf('[跳过统计] %s - %s (seed=%d): 文件不存在\n', dname, agg, s);
                continue;
            end

            d_dac = load(dac_file);
            d_ac  = load(ac_file);

            iters_dac_all(end+1) = d_dac.iter_converge;
            comm_dac_all(end+1)  = d_dac.comm_train;
            final_dac_all(end+1) = d_dac.conv_curve_dac(end);

            iters_ac_all(end+1)  = d_ac.iter_converge;
            comm_ac_all(end+1)   = d_ac.comm_train;
            final_ac_all(end+1)  = d_ac.conv_curve_ac(end);
        end

        if isempty(iters_dac_all) || isempty(iters_ac_all)
            continue;
        end

        % 均值 ± 标准差（标准差为0时, 即只有1个seed, std返回0）
        iters_dac_mean = mean(iters_dac_all); iters_dac_std = std(iters_dac_all);
        comm_dac_mean  = mean(comm_dac_all);  comm_dac_std  = std(comm_dac_all);
        final_dac_mean = mean(final_dac_all); final_dac_std = std(final_dac_all);

        iters_ac_mean = mean(iters_ac_all); iters_ac_std = std(iters_ac_all);
        comm_ac_mean  = mean(comm_ac_all);  comm_ac_std  = std(comm_ac_all);
        final_ac_mean = mean(final_ac_all); final_ac_std = std(final_ac_all);

        summary_table(end+1,:) = {dname, agg_labels{a_idx}, ...
            iters_dac_mean, iters_dac_std, comm_dac_mean, comm_dac_std, final_dac_mean, final_dac_std, ...
            iters_ac_mean,  iters_ac_std,  comm_ac_mean,  comm_ac_std,  final_ac_mean,  final_ac_std, ...
            numel(iters_dac_all)};
    end

    xlabel('Iteration'); ylabel('Max agent disagreement');
    title(dname);
    legend(legend_entries, 'Location','northeast', 'FontSize',7, 'NumColumns',2);
    grid on;
    hold off;
end

sgtitle(sprintf('Consensus Convergence Summary (Train=%.0f%%, curves: seed=%d, stats: seed=%d~%d)', ...
    train_ratio*100, plot_seed, seed_list(1), seed_list(end)));

%% 命令行汇总表（均值 ± 标准差，跨 seed_list）
fprintf('\n==========================================================================================================================\n');
fprintf('  Consensus Convergence & ET Communication Summary  (Train=%.0f%%, averaged over %d seeds: %d~%d)\n', ...
    train_ratio*100, numel(seed_list), seed_list(1), seed_list(end));
fprintf('==========================================================================================================================\n');
fprintf('%-12s %-6s %6s | %14s %14s %14s | %14s %14s %14s\n', ...
    'Dataset','Agg','N', 'DAC_it', 'DAC_ET', 'DAC_final', 'AC_it', 'AC_ET', 'AC_final');
fprintf('%s\n', repmat('-',1,122));
for i = 1:size(summary_table,1)
    row = summary_table(i,:);
    fprintf('%-12s %-6s %6d | %6.1f%+6.1f  %6.1f%+6.1f  %6.2e%+6.1e | %6.1f%+6.1f  %6.1f%+6.1f  %6.2e%+6.1e\n', ...
        row{1}, row{2}, row{15}, ...
        row{3}, row{4}, row{5}, row{6}, row{7}, row{8}, ...
        row{9}, row{10}, row{11}, row{12}, row{13}, row{14});
end
fprintf('%s\n\n', repmat('-',1,122));

%% 保存汇总表（方便后续直接读取，不用重新跑）
save(fullfile('Result','Dataset','consensus_summary_allseeds.mat'), 'summary_table', 'seed_list', 'plot_seed', 'datasets', 'agg_list', 'agg_labels', 'train_ratio');
fprintf('汇总表已保存到 Result/Dataset/consensus_summary_allseeds.mat\n');

%% 导出 CSV，方便直接粘贴进 Word/Excel 表格
csv_path = fullfile('Result','Dataset','consensus_summary_allseeds.csv');
fid = fopen(csv_path, 'w');
fprintf(fid, 'Dataset,Agg,N,DAC_iters_mean,DAC_iters_std,DAC_ET_mean,DAC_ET_std,DAC_final_mean,DAC_final_std,AC_iters_mean,AC_iters_std,AC_ET_mean,AC_ET_std,AC_final_mean,AC_final_std\n');
for i = 1:size(summary_table,1)
    row = summary_table(i,:);
    fprintf(fid, '%s,%s,%d,%.4f,%.4f,%.4f,%.4f,%.6e,%.6e,%.4f,%.4f,%.4f,%.4f,%.6e,%.6e\n', ...
        row{1}, row{2}, row{15}, row{3}, row{4}, row{5}, row{6}, row{7}, row{8}, ...
        row{9}, row{10}, row{11}, row{12}, row{13}, row{14});
end
fclose(fid);

%% Final ET consensus convergence values
FigureObj_FinalConsensus = figure('Color','w','Position',[80,80,1300,850], ...
    'Name','Final ET consensus convergence values');

for d_idx = 1:numel(datasets)
    dname = datasets{d_idx};
    ax = subplot(2,2,d_idx,'Parent',FigureObj_FinalConsensus);
    hold(ax,'on');

    final_dac = nan(1,numel(agg_labels));
    final_ac  = nan(1,numel(agg_labels));
    for a_idx = 1:numel(agg_labels)
        row_idx = find(strcmp(summary_table(:,1), dname) & ...
            strcmp(summary_table(:,2), agg_labels{a_idx}), 1);
        if ~isempty(row_idx)
            final_dac(a_idx) = summary_table{row_idx,7};
            final_ac(a_idx)  = summary_table{row_idx,13};
        end
    end

    bar(ax, [final_dac(:), final_ac(:)]);
    set(ax, 'YScale', 'log');
    xticks(ax, 1:numel(agg_labels));
    xticklabels(ax, agg_labels);
    ylabel(ax, 'Final consensus disagreement');
    title(ax, dname);
    grid(ax,'on');
    legend(ax, {'IP-DAC','IP-AC'}, 'Location','best');
end

sgtitle(FigureObj_FinalConsensus, ...
    sprintf('Final ET consensus convergence values (Train=%.0f%%, seeds=%d-%d)', ...
    train_ratio*100, seed_list(1), seed_list(end)));

final_fig_path = fullfile('Result','Dataset','consensus_final_values.fig');
final_png_path = fullfile('Result','Dataset','consensus_final_values.png');
savefig(FigureObj_FinalConsensus, final_fig_path);
saveas(FigureObj_FinalConsensus, final_png_path);
fprintf('Final convergence figure saved to %s and %s\n', final_fig_path, final_png_path);
fprintf('CSV 已导出到 %s，可直接用 Excel 打开或粘贴进 Word 表格\n', csv_path);
