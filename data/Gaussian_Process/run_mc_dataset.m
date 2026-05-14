% run_mc_dataset.m

clc; clear; close all;


diary('Overnight_Log.txt');
diary on;
fprintf('\n======================================================\n');
fprintf('  任务开始时间: %s\n', datestr(now));
fprintf('======================================================\n');

%%  1. 全局配置区域 =================
datasets    ={'KIN40K','POL','PUMADYN32NM','SARCOS'};
%datasets    ={'KIN40K'};
train_ratio = 0.4;
n_mc        = 5;   

run_log_and_ip = true;  % 跑 LoG-GP基线 和 诱导点 (DAC+AC)
run_tp         = true;  % 跑 测试点 (DAC+AC)
run_nbr        = true;  % 跑 邻域静态聚合

%% ================= 2. 方法 =================
methods_dict = {
    'LoG-MOE',      'log_moe_mc%d.mat';           'LoG-GPOE',     'log_gpoe_mc%d.mat';
    'IP-DAC-MOE',   'moe_tr%d_mc%d.mat';          'IP-DAC-GPOE',  'gpoe_tr%d_mc%d.mat';
    'IP-DAC-POE',   'poe_tr%d_mc%d.mat';          'IP-DAC-BCM',   'bcm_tr%d_mc%d.mat';
    'IP-DAC-RBCM',  'rbcm_tr%d_mc%d.mat';         'IP-AC-MOE',    'moe_ac_tr%d_mc%d.mat';
    'IP-AC-GPOE',   'gpoe_ac_tr%d_mc%d.mat';      'IP-AC-POE',    'poe_ac_tr%d_mc%d.mat';
    'IP-AC-BCM',    'bcm_ac_tr%d_mc%d.mat';       'IP-AC-RBCM',   'rbcm_ac_tr%d_mc%d.mat';
    'TP-DAC-MOE',   'moe_tp_tr%d_mc%d.mat';       'TP-DAC-GPOE',  'gpoe_tp_tr%d_mc%d.mat';
    'TP-DAC-POE',   'poe_tp_tr%d_mc%d.mat';       'TP-DAC-BCM',   'bcm_tp_tr%d_mc%d.mat';
    'TP-DAC-RBCM',  'rbcm_tp_tr%d_mc%d.mat';      'TP-AC-MOE',    'moe_ac_tp_tr%d_mc%d.mat';
    'TP-AC-GPOE',   'gpoe_ac_tp_tr%d_mc%d.mat';   'TP-AC-POE',    'poe_ac_tp_tr%d_mc%d.mat';
    'TP-AC-BCM',    'bcm_ac_tp_tr%d_mc%d.mat';    'TP-AC-RBCM',   'rbcm_ac_tp_tr%d_mc%d.mat';
    'NBR-MOE',      'moe_nbr_mc%d.mat';      'NBR-GPOE',     'gpoe_nbr_mc%d.mat';
    'NBR-POE',      'poe_nbr_mc%d.mat';      'NBR-BCM',      'bcm_nbr_mc%d.mat';
    'NBR-RBCM',     'rbcm_nbr_mc%d.mat';
};
n_methods = size(methods_dict, 1);
tr_tag = round(train_ratio * 100);

%% ================= 3. 执行仿真 (做菜阶段) =================
for d = 1:numel(datasets)
    for mc = 1:n_mc
        fprintf('\n\n========== 当前进度: %s (MC %d/%d) ==========\n', datasets{d}, mc, n_mc);
        
        % 任务 1：LoG 和 诱导点
        if run_log_and_ip
            try
                run_LoGGP_comparison(datasets{d}, train_ratio, mc);
                run_inducingpoint_dataset(datasets{d}, 'all', train_ratio, mc);
            catch ME, fprintf('\n[错误跳过] 诱导点任务失败: %s\n', ME.message); end
        end
        
        % 任务 2：测试点
        if run_tp
            try
                run_testpoint_dataset(datasets{d}, 'all', train_ratio, mc);
            catch ME, fprintf('\n[错误跳过] 测试点任务失败: %s\n', ME.message); end
        end
        
        % 任务 3：邻域静态聚合
        if run_nbr
            try
                run_neighbor_dataset(datasets{d}, 'all', train_ratio, mc);
            catch ME, fprintf('\n[错误跳过] 邻域任务失败: %s\n', ME.message); end
        end
        
        % 每次循环结束，强制清理一次内存，保证 5GB 内存绝对不爆！
        clearvars -except datasets train_ratio n_mc run_log_and_ip run_tp run_nbr methods_dict n_methods tr_tag d mc;
        pause(1); 
    end
end

%% ================= 4. 读取与汇总 (端菜上桌阶段) =================
fprintf('\n\n======================================================\n');
fprintf('  仿真阶段结束，开始从硬盘读取 .mat 文件并生成汇总表...\n');
fprintf('======================================================\n');

mean_results = nan(numel(datasets), n_methods, 5); % SMSE, RMSE, NLPD, TrT, TeT
std_results  = nan(numel(datasets), n_methods, 5);

for d = 1:numel(datasets)
    SaveFolder = fullfile('Result', 'Dataset', datasets{d});
    for mi = 1:n_methods
        file_template = methods_dict{mi, 2};
        raw = nan(n_mc, 5);
        for mc = 1:n_mc
            if contains(file_template, 'log_') || contains(file_template, '_nbr_')
                fname = fullfile(SaveFolder, sprintf(file_template, mc));
            else
                fname = fullfile(SaveFolder, sprintf(file_template, tr_tag, mc));
            end
            if exist(fname, 'file')
                r = load(fname, 'smse', 'rmse', 'nlpd', 't_train', 't_test');
                raw(mc, :) = [r.smse, r.rmse, r.nlpd, r.t_train, r.t_test];
            end
        end
        mean_results(d, mi, :) = mean(raw, 1, 'omitnan');
        std_results(d, mi, :)  = std(raw, 0, 1, 'omitnan');
    end
end

%% ================= 5. 打印超长精简大表 =================
group_idx = [1, 3, 8, 13, 18, 23, 28]; % 分组边界

for d = 1:numel(datasets)
    fprintf('\n\n==================== %s (Train=%.0f%%, MC=%d) ====================\n', datasets{d}, train_ratio*100, n_mc);
    fprintf('  %-15s  %11s  %11s  %11s  %10s  %10s\n', 'Method', 'SMSE', 'RMSE', 'NLPD', 'Train_T(s)', 'Test_T(s)');
    fprintf('  %s\n', repmat('-', 1, 80));
    
    for g = 1:6
        group_means = mean_results(d, group_idx(g):group_idx(g+1)-1, 1);
        if all(isnan(group_means(:))), continue; end % 如果这一组没跑，跳过打印
        
        for mi = group_idx(g) : group_idx(g+1)-1
            m = squeeze(mean_results(d, mi, :));
            s = squeeze(std_results(d, mi, :));
            if all(isnan(m)), continue; end
            
            fprintf('  %-15s  %.4f±%.4f  %.4f±%.4f  %5.2f±%4.2f  %7.2f±%4.2f  %7.3f±%4.3f\n', ...
                methods_dict{mi, 1}, m(1), s(1), m(2), s(2), m(3), s(3), m(4), s(4), m(5), s(5));
        end
        if g < 6, fprintf('  %s\n', repmat('-', 1, 80)); end 
    end
    fprintf('=================================================================================\n');
end

% 保存汇总矩阵
save(fullfile('Result', 'Dataset', 'All_27_Methods_Summary.mat'), 'datasets', 'methods_dict', 'mean_results', 'std_results');
fprintf('\n汇总完成，结果矩阵已保存至 All_27_Methods_Summary.mat\n');

fprintf('\n夜间挂机任务结束时间: %s\n', datestr(now));
diary off;
