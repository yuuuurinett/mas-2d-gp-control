% run_mc_dataset.m
clc; clear; close all;

diary('Overnight_Log.txt');
diary on;
fprintf('\n======================================================\n');
fprintf('  任务开始时间: %s\n', datestr(now));
fprintf('======================================================\n');

%%  1. 全局配置区域
datasets    = {'KIN40K','POL','PUMADYN32NM','SARCOS'};
%datasets    = {'SARCOS'};
train_ratio = 0.4;
n_mc        = 3;

run_log        = true;
run_ip         = true;
run_tp         = true;
run_cen        = true;
run_nbr        = true;

%% 2. 方法字典
methods_dict = {
    'LoG-MOE',      'log_moe_mc%d.mat';           'LoG-GPOE',     'log_gpoe_mc%d.mat';
    'IP-DAC-MOE',   'moe_tr%d_mc%d.mat';          'IP-DAC-GPOE',  'gpoe_tr%d_mc%d.mat';
    'IP-DAC-POE',   'poe_tr%d_mc%d.mat';          'IP-DAC-BCM',   'bcm_tr%d_mc%d.mat';
    'IP-DAC-RBCM',  'rbcm_tr%d_mc%d.mat';
    'IP-AC-MOE',    'moe_ac_tr%d_mc%d.mat';       'IP-AC-GPOE',   'gpoe_ac_tr%d_mc%d.mat';
    'IP-AC-POE',    'poe_ac_tr%d_mc%d.mat';       'IP-AC-BCM',    'bcm_ac_tr%d_mc%d.mat';
    'IP-AC-RBCM',   'rbcm_ac_tr%d_mc%d.mat';
    'TP-DAC-MOE',   'moe_tp_tr%d_mc%d.mat';       'TP-DAC-GPOE',  'gpoe_tp_tr%d_mc%d.mat';
    'TP-DAC-POE',   'poe_tp_tr%d_mc%d.mat';       'TP-DAC-BCM',   'bcm_tp_tr%d_mc%d.mat';
    'TP-DAC-RBCM',  'rbcm_tp_tr%d_mc%d.mat';
    'TP-AC-MOE',    'moe_ac_tp_tr%d_mc%d.mat';    'TP-AC-GPOE',   'gpoe_ac_tp_tr%d_mc%d.mat';
    'TP-AC-POE',    'poe_ac_tp_tr%d_mc%d.mat';    'TP-AC-BCM',    'bcm_ac_tp_tr%d_mc%d.mat';
    'TP-AC-RBCM',   'rbcm_ac_tp_tr%d_mc%d.mat';
    'CEN-MOE',      'moe_cen_tr%d_mc%d.mat';      'CEN-GPOE',     'gpoe_cen_tr%d_mc%d.mat';
    'CEN-POE',      'poe_cen_tr%d_mc%d.mat';      'CEN-BCM',      'bcm_cen_tr%d_mc%d.mat';
    'CEN-RBCM',     'rbcm_cen_tr%d_mc%d.mat';
    'NBR-MOE',      'moe_nbr_tr%d_mc%d.mat';      'NBR-GPOE',     'gpoe_nbr_tr%d_mc%d.mat';
    'NBR-POE',      'poe_nbr_tr%d_mc%d.mat';      'NBR-BCM',      'bcm_nbr_tr%d_mc%d.mat';
    'NBR-RBCM',     'rbcm_nbr_tr%d_mc%d.mat';
};
methods_names = methods_dict(:,1);
methods_files = methods_dict(:,2);
num_methods   = length(methods_names);

% 8 种指标: SMSE, RMSE, Train(ms/pt), Test(ms/pt), Comm_Tr, Comm_Te, Iter_Conv
mean_results = NaN(length(datasets), num_methods, 7);
std_results  = NaN(length(datasets), num_methods, 7);
tr_tag       = round(train_ratio * 100);

%% 3. 运行 & 统计
for d = 1:length(datasets)
    dname = datasets{d};
    fprintf('\n>>> 处理数据集: %s <<<\n', dname);

    for seed = 1:n_mc
        if run_log, run_LoGGP_comparison(dname, train_ratio, seed); end
        if run_ip,  run_inducingpoint_dataset(dname, 'all', train_ratio, seed); end
        if run_tp,  run_testpoint_dataset(dname, 'all', train_ratio, seed); end
        if run_cen, run_centralized_dataset(dname, 'all', train_ratio, seed); end
        if run_nbr, run_neighbor_dataset(dname, 'all', train_ratio, seed); end
    end

    for mi = 1:num_methods
        sm=NaN(1,n_mc); rm=NaN(1,n_mc); %nl=NaN(1,n_mc);
        t_tr=NaN(1,n_mc); t_te=NaN(1,n_mc);
        c_tr=NaN(1,n_mc); c_te=NaN(1,n_mc); it=NaN(1,n_mc);

        for mc = 1:n_mc
            if contains(methods_files{mi}, 'tr%d')
                fname = sprintf(methods_files{mi}, tr_tag, mc);
            else
                fname = sprintf(methods_files{mi}, mc);
            end
            file_path = fullfile('Result', 'Dataset', dname, fname);

            if exist(file_path, 'file')
                try
                    res = load(file_path);
                    sm(mc)   = res.smse;
                    rm(mc)   = res.rmse;
                    %nl(mc)   = res.nlpd;
                    t_tr(mc) = res.t_train_per_point;
                    t_te(mc) = res.t_test_per_point;
                    if isfield(res,'comm_train'),   c_tr(mc) = res.comm_train;   else, c_tr(mc) = 0; end
                    if isfield(res,'comm_test'),    c_te(mc) = res.comm_test;    else, c_te(mc) = 0; end
                    if isfield(res,'iter_converge'), it(mc)  = res.iter_converge; else, it(mc)  = NaN; end
                catch
                    fprintf('  [读取失败] %s\n', fname);
                end
            end
        end
        mean_results(d,mi,1)=mean(sm,'omitnan');
        std_results(d,mi,1)=std(sm,'omitnan');

        mean_results(d,mi,2)=mean(rm,'omitnan');
        std_results(d,mi,2)=std(rm,'omitnan');

        mean_results(d,mi,3)=mean(t_tr,'omitnan');
        std_results(d,mi,3)=std(t_tr,'omitnan');

        mean_results(d,mi,4)=mean(t_te,'omitnan');
        std_results(d,mi,4)=std(t_te,'omitnan');

        mean_results(d,mi,5)=mean(c_tr,'omitnan');
        std_results(d,mi,5)=0;

        mean_results(d,mi,6)=mean(c_te,'omitnan');
        std_results(d,mi,6)=0;

        mean_results(d,mi,7)=mean(it,'omitnan');
        std_results(d,mi,7)=0;
    end
end

%% 4. 打印排版
metrics_names = {'SMSE','RMSE','Train_T(ms/pt)','Test_T(ms/pt)',...
                 'Comm_Train','Comm_Test','Iter_Converge'};
agg_list     = {'MOE','GPOE','POE','BCM','RBCM'};
col_prefixes = {'LoG','CEN','IP-DAC','IP-AC','TP-DAC','TP-AC','NBR'};
sep_wide = repmat('=',1,130);
sep_thin = repmat('-',1,130);

for met = 1:7
    fprintf('\n%s\n', sep_wide);
    fprintf('  Metric: %s  (Train=%.0f%%  MC=%d)\n', metrics_names{met}, train_ratio*100, n_mc);
    fprintf('%s\n', sep_wide);
    fprintf('  Agg    Dataset                   LoG             CEN          IP-DAC           IP-AC          TP-DAC           TP-AC             NBR\n');
    fprintf('  %s\n', sep_thin);

    for a_idx = 1:length(agg_list)
        agg = agg_list{a_idx};
        for d = 1:length(datasets)
            if d==1
                fprintf('  %-4s   %-12s', agg, datasets{d});
            else
                fprintf('         %-12s', datasets{d});
            end

            for col = 1:7
                prefix      = col_prefixes{col};
                target_name = sprintf('%s-%s', prefix, agg);
                mi = find(strcmp(methods_names, target_name));

                if isempty(mi)
                    fprintf('  %14s', '-');
                else
                    mv = mean_results(d,mi,met);
                    sv = std_results(d,mi,met);
                    if isnan(mv)
                        fprintf('  %14s', '-');
                    elseif met >= 5
                        % Comm 和 Iter 列：整数显示，0 显示为 -
                        if mv==0 || isnan(mv)
                            fprintf('  %14s', '-');
                        else
                            fprintf('  %14.0f', mv);
                        end
                    elseif met==3 || met==4 
                        fprintf('  %14s', sprintf('%.2f±%.2f', mv, sv));
                    else
                        fprintf('  %14s', sprintf('%.4f±%.4f', mv, sv));
                    end
                end
            end
            fprintf('\n');
        end
        fprintf('  %s\n', sep_thin);
    end
end

fprintf('\n%s\n', sep_wide);
fprintf('  任务结束时间: %s\n', datestr(now));
fprintf('%s\n', sep_wide);

save(fullfile('Result','Dataset','All_Methods_Summary.mat'), 'mean_results', 'std_results');
diary off;