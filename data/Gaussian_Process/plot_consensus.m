%% plot_consensus_convergence.m
% 读取mat文件，画：
% 1. IP-DAC 和 IP-AC 的收敛曲线（disagreement vs iter）
% 2. ET触发次数 vs 总迭代步数（通信节省对比）
clc; close all;

DatasetName = 'KIN40K';
train_ratio = 0.4;
seed        = 1;
tr_tag      = round(train_ratio*100);
agg_method  = 'poe';
SaveFolder  = fullfile('Result', 'Dataset', DatasetName);

%% 读取mat文件
dac_file = fullfile(SaveFolder, sprintf('%s_tr%d_mc%d.mat',    agg_method, tr_tag, seed));
ac_file  = fullfile(SaveFolder, sprintf('%s_ac_tr%d_mc%d.mat', agg_method, tr_tag, seed));

if ~exist(dac_file,'file') || ~exist(ac_file,'file')
    error('找不到mat文件，请先运行 run_inducingpoint_dataset');
end

d_dac = load(dac_file);
d_ac  = load(ac_file);

%% 命令行输出详细数值
fprintf('\n========== Consensus Convergence Summary ==========\n');
fprintf('Dataset: %s  Method: %s  seed=%d\n', DatasetName, upper(agg_method), seed);
fprintf('\n');

iters_dac  = d_dac.iter_converge;
iters_ac   = d_ac.iter_converge;
comm_dac   = d_dac.comm_train;
comm_ac    = d_ac.comm_train;
final_dac  = d_dac.conv_curve_dac(end);
final_ac   = d_ac.conv_curve_ac(end);
saving_dac = (1 - comm_dac / iters_dac) * 100;
saving_ac  = (1 - comm_ac  / iters_ac)  * 100;

fprintf('%-24s  %10s  %10s\n', '', 'IP-DAC', 'IP-AC');
fprintf('%s\n', repmat('-',1,48));
fprintf('%-24s  %10d  %10d\n',   '总迭代步数',         iters_dac, iters_ac);
fprintf('%-24s  %10.1f  %10.1f\n','ET触发次数/agent',   comm_dac,  comm_ac);
fprintf('%-24s  %9.1f%%  %9.1f%%\n','通信节省比例',     saving_dac, saving_ac);
fprintf('%-24s  %10.2e  %10.2e\n','最终disagreement(DAC)', final_dac, final_ac);
fprintf('%s\n\n', repmat('-',1,48));

%% 画图
figure('Color','w','Position',[100,100,1100,430]);

%% 左图：收敛曲线，分别标注各自最终收敛值
subplot(1,2,1);
semilogy(1:iters_dac, d_dac.conv_curve_dac, 'b-',  'LineWidth',2.0); hold on;
semilogy(1:iters_ac,  d_ac.conv_curve_ac,   'r--', 'LineWidth',2.0);

% 分别用虚线标注各自最终收敛值
yline(final_dac, 'b:', 'LineWidth',1.5);
yline(final_ac,  'r:', 'LineWidth',1.5);

% 右侧标注具体数值
xlim_max = max(iters_dac, iters_ac);
text(xlim_max*0.98, final_dac*2, sprintf('DAC: %.2e', final_dac), ...
    'Color','b','FontSize',9,'FontWeight','bold','HorizontalAlignment','right');
text(xlim_max*0.98, final_ac*0.5, sprintf('AC: %.2e', final_ac), ...
    'Color','r','FontSize',9,'FontWeight','bold','HorizontalAlignment','right');

xlabel('Iteration'); ylabel('Max agent disagreement');
title(sprintf('Consensus Convergence [%s, %s]', DatasetName, upper(agg_method)));
legend(sprintf('IP-DAC (final=%.2e)', final_dac), ...
       sprintf('IP-AC  (final=%.2e)', final_ac), ...
       'Location','northeast');
grid on;

%% 右图：ET触发次数 vs 总迭代步数
subplot(1,2,2);
data = [iters_dac, comm_dac; iters_ac, comm_ac];
b = bar(data, 'grouped');
b(1).FaceColor = [0.7 0.7 0.7];
b(2).FaceColor = [0.2 0.5 0.9];

for k = 1:2
    for g = 1:2
        text(k + (g-1.5)*0.22, data(k,g) + max(data(:))*0.01, ...
            sprintf('%.0f', data(k,g)), ...
            'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
    end
end

text(1, max(data(1,:))*1.12, sprintf('节省 %.1f%%', saving_dac), ...
    'HorizontalAlignment','center','FontSize',10,'Color',[0 0.5 0],'FontWeight','bold');
text(2, max(data(2,:))*1.12, sprintf('节省 %.1f%%', saving_ac), ...
    'HorizontalAlignment','center','FontSize',10,'Color',[0 0.5 0],'FontWeight','bold');

set(gca,'XTickLabel',{'IP-DAC','IP-AC'});
legend('总迭代步数','ET触发次数','Location','northeast');
ylabel('次数');
title('ET Communication Savings');
ylim([0, max(data(:))*1.25]);
grid on;