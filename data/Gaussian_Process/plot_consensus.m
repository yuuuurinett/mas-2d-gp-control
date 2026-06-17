%% plot_consensus_convergence_summary.m
% 汇总图：4个数据集 x 5个方法的收敛曲线对比
% 每个数据集一个子图，每个方法用不同颜色展示 IP-DAC（实线）和 IP-AC（虚线）
clc; close all;

datasets    = {'KIN40K','POL','PUMADYN32NM','SARCOS'};
agg_list    = {'moe','gpoe','poe','bcm','rbcm'};
agg_labels  = {'MOE','GPOE','POE','BCM','RBCM'};
train_ratio = 0.4;
seed        = 1;
tr_tag      = round(train_ratio*100);

colors = lines(5);  % 5种颜色对应5个method

figure('Color','w','Position',[60,60,1400,950]);

summary_table = {};  % 汇总数值表

for d_idx = 1:length(datasets)
    dname = datasets{d_idx};
    SaveFolder = fullfile('Result','Dataset',dname);

    subplot(2,2,d_idx);
    hold on;
    legend_entries = {};

    for a_idx = 1:length(agg_list)
        agg = agg_list{a_idx};
        dac_file = fullfile(SaveFolder, sprintf('%s_tr%d_mc%d.mat',    agg, tr_tag, seed));
        ac_file  = fullfile(SaveFolder, sprintf('%s_ac_tr%d_mc%d.mat', agg, tr_tag, seed));

        if ~exist(dac_file,'file') || ~exist(ac_file,'file')
            fprintf('[跳过] %s - %s: 文件不存在\n', dname, agg);
            continue;
        end

        d_dac = load(dac_file);
        d_ac  = load(ac_file);

        iters_dac = d_dac.iter_converge;
        iters_ac  = d_ac.iter_converge;
        final_dac = d_dac.conv_curve_dac(end);
        final_ac  = d_ac.conv_curve_ac(end);
        comm_dac  = d_dac.comm_train;
        comm_ac   = d_ac.comm_train;

        semilogy(1:iters_dac, d_dac.conv_curve_dac, '-',  'Color', colors(a_idx,:), 'LineWidth',1.6);
        semilogy(1:iters_ac,  d_ac.conv_curve_ac,   '--', 'Color', colors(a_idx,:), 'LineWidth',1.6);

        legend_entries{end+1} = sprintf('%s-DAC', agg_labels{a_idx});
        legend_entries{end+1} = sprintf('%s-AC',  agg_labels{a_idx});

        summary_table(end+1,:) = {dname, agg_labels{a_idx}, iters_dac, comm_dac, final_dac, ...
                                          iters_ac,  comm_ac,  final_ac};
    end

    xlabel('Iteration'); ylabel('Max agent disagreement');
    title(dname);
    legend(legend_entries, 'Location','northeast', 'FontSize',7, 'NumColumns',2);
    grid on;
    hold off;
end

sgtitle(sprintf('Consensus Convergence Summary (Train=%.0f%%, seed=%d)', train_ratio*100, seed));

%% 命令行汇总表
fprintf('\n================================================================================================\n');
fprintf('  Consensus Convergence & ET Communication Summary  (Train=%.0f%%, seed=%d)\n', train_ratio*100, seed);
fprintf('================================================================================================\n');
fprintf('%-12s %-6s %8s %8s %10s | %8s %8s %10s\n', ...
    'Dataset','Agg','DAC_it','DAC_ET','DAC_final','AC_it','AC_ET','AC_final');
fprintf('%s\n', repmat('-',1,96));
for i = 1:size(summary_table,1)
    row = summary_table(i,:);
    fprintf('%-12s %-6s %8d %8.1f %10.2e | %8d %8.1f %10.2e\n', row{:});
end
fprintf('%s\n\n', repmat('-',1,96));