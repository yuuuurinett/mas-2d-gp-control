%% plot_consensus_convergence.m
% 直接从保存的mat文件读取收敛曲线，无需重跑仿真
clc; close all;

DatasetName = 'KIN40K';
train_ratio = 0.4;
seed = 1;
tr_tag = round(train_ratio*100);
SaveFolder = fullfile('Result','Dataset',DatasetName);

% 读取DAC方法（如poe）
dac_file = fullfile(SaveFolder, sprintf('poe_tr%d_mc%d.mat', tr_tag, seed));
ac_file  = fullfile(SaveFolder, sprintf('poe_ac_tr%d_mc%d.mat', tr_tag, seed));

if ~exist(dac_file,'file') || ~exist(ac_file,'file')
    error('找不到mat文件，请先运行 run_inducingpoint_dataset');
end

d_dac = load(dac_file, 'conv_curve_dac', 'iter_converge', 'comm_train');
d_ac  = load(ac_file,  'conv_curve_ac',  'iter_converge', 'comm_train');

figure('Color','w','Position',[100,100,700,350]);
semilogy(1:length(d_dac.conv_curve_dac), d_dac.conv_curve_dac, 'b-',  'LineWidth',1.5); hold on;
semilogy(1:length(d_ac.conv_curve_ac),   d_ac.conv_curve_ac,   'r--', 'LineWidth',1.5);
xlabel('Iteration'); ylabel('Max agent disagreement');
title(sprintf('Consensus Convergence [%s, POE]', DatasetName));
legend(sprintf('IP-DAC (ET=%.0f)', d_dac.comm_train), ...
       sprintf('IP-AC  (ET=%.0f)', d_ac.comm_train), ...
       'Location','northeast');
grid on;