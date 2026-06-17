% run_mc_dataset.m
clc; clear; close all;

diary('Overnight_Log.txt');
diary on;
fprintf('\n======================================================\n');
fprintf('  任务开始时间: %s\n', datestr(now));
fprintf('======================================================\n');

%%  1. 全局配置区域
datasets    = {'KIN40K','POL','PUMADYN32NM','SARCOS'};
train_ratio = 0.4;
n_mc        = 1;

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
        sm=NaN(1,n_mc); rm=NaN(1,n_mc);
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

%% 4. 🌟 升级版排版区域：同时打印控制台并生成完美 Word HTML 大表
metrics_names = {'SMSE','RMSE','Train_T(ms/pt)','Test_T(ms/pt)',...
                 'Comm_Train','Comm_Test','Iter_Converge'};
agg_list     = {'MOE','GPOE','POE','BCM','RBCM'};
col_prefixes = {'LoG','CEN','IP-DAC','IP-AC','TP-DAC','TP-AC','NBR'};
sep_wide = repmat('=',1,130);
sep_thin = repmat('-',1,130);

% 打开 HTML 文件流
html_file = 'All_Metrics_Tables_For_Word.html';
fid_html = fopen(html_file, 'w', 'n', 'utf-8');

% 写入高级、紧凑、防换行的样式配置
fprintf(fid_html, '<!DOCTYPE html><html><head><style>\n');
fprintf(fid_html, 'h3 { font-family: Arial, sans-serif; margin-top: 30px; color: #333; }\n');
fprintf(fid_html, 'table { border-collapse: collapse; font-family: "Times New Roman", Times, serif; font-size: 10.5pt; text-align: center; margin-bottom: 20px; }\n');
fprintf(fid_html, 'th, td { border: 1px solid black; padding: 4px 6px; white-space: nowrap; }\n');
fprintf(fid_html, 'th { background-color: #DDEBF7; font-weight: bold; }\n'); % 经典淡蓝色表头
fprintf(fid_html, '</style></head><body>\n');

for met = 1:7
    % --- 控制台标准打印 ---
    fprintf('\n%s\n', sep_wide);
    fprintf('  Metric: %s  (Train=%.0f%%  MC=%d)\n', metrics_names{met}, train_ratio*100, n_mc);
    fprintf('%s\n', sep_wide);
    fprintf('  Agg    Dataset                   LoG             CEN          IP-DAC           IP-AC          TP-DAC           TP-AC             NBR\n');
    fprintf('  %s\n', sep_thin);

    % --- HTML 表头写入 ---
    fprintf(fid_html, '<h3>Metric: %s &nbsp;(Train=%.0f%%, MC=%d)</h3>\n', metrics_names{met}, train_ratio*100, n_mc);
    fprintf(fid_html, '<table>\n');
    fprintf(fid_html, '  <tr><th>聚合</th><th>数据集</th><th>LoG</th><th>CEN</th><th>IP-DAC</th><th>IP-AC</th><th>TP-DAC</th><th>TP-AC</th><th>NBR</th></tr>\n');

    for a_idx = 1:length(agg_list)
        agg = agg_list{a_idx};
        for d = 1:length(datasets)
            % 控制台每行开头处理
            if d==1, fprintf('  %-4s   %-12s', agg, datasets{d});
            else,    fprintf('         %-12s', datasets{d}); end

            % HTML 每行开头处理
            fprintf(fid_html, '  <tr>\n');
            if d == 1
                % 🌟 核心：第一行自动纵向合并单元格，完美复刻模板效果
                fprintf(fid_html, '    <td rowspan="%d" style="vertical-align: middle; font-weight: bold; background-color: #F2F2F2;">%s</td>\n', length(datasets), agg);
            end
            fprintf(fid_html, '    <td style="font-weight: bold;">%s</td>\n', datasets{d});

            % 遍历 7 个方法列
            for col = 1:7
                prefix      = col_prefixes{col};
                target_name = sprintf('%s-%s', prefix, agg);
                mi = find(strcmp(methods_names, target_name));

                if isempty(mi)
                    fprintf('  %14s', '-');
                    fprintf(fid_html, '    <td>-</td>\n');
                else
                    mv = mean_results(d,mi,met);
                    sv = std_results(d,mi,met);
                    
                    if isnan(mv)
                        fprintf('  %14s', '-');
                        fprintf(fid_html, '    <td>-</td>\n');
                    elseif met >= 5
                        if mv==0
                            fprintf('  %14s', '-');
                            fprintf(fid_html, '    <td>-</td>\n');
                        else
                            fprintf('  %14.0f', mv);
                            fprintf(fid_html, '    <td>%.0f</td>\n', mv);
                        end
                    elseif met==3 || met==4 
                        str_val = sprintf('%.2f±%.2f', mv, sv);
                        fprintf('  %14s', str_val);
                        fprintf(fid_html, '    <td>%s</td>\n', str_val); % 🌟 掐掉空格，杜绝换行
                    else
                        str_val = sprintf('%.4f±%.4f', mv, sv);
                        fprintf('  %14s', str_val);
                        fprintf(fid_html, '    <td>%s</td>\n', str_val); % 🌟 掐掉空格，杜绝换行
                    end
                end
            end
            fprintf('\n');
            fprintf(fid_html, '  </tr>\n');
        end
        fprintf('  %s\n', sep_thin);
    end
    fprintf(fid_html, '</table>\n');
end

fprintf('\n%s\n', sep_wide);
fprintf('  任务结束时间: %s\n', datestr(now));
fprintf('%s\n', sep_wide);

fprintf(fid_html, '</body></html>\n');
fclose(fid_html);

save(fullfile('Result','Dataset','All_Methods_Summary.mat'), 'mean_results', 'std_results');
diary off;

fprintf('\n📊 完美 Word 表格已同步生成！\n');
fprintf('请在当前文件夹双击打开【All_Metrics_Tables_For_Word.html】\n');
fprintf('然后：Ctrl+A (全选) -> Ctrl+C (复制) -> 直接粘贴进你的 Word 文档！\n\n');