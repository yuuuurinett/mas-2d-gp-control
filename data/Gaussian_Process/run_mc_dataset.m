%% run_mc_dataset_markdown_transposed.m
% Robust MC runner + transposed Markdown report.
%
% Main purpose:
%   The original command-window table is too wide for Word.
%   This version additionally exports a Markdown report with the axes swapped:
%
%       For each metric:
%           For each aggregation method MOE/GPOE/POE/BCM/RBCM:
%               rows    = LoG, CEN, IP-DAC, IP-AC, TP-DAC, TP-AC, NBR
%               columns = datasets
%
% Output:
%   Result/Dataset/All_Methods_Summary.mat
%   Result/Dataset/All_Methods_Summary_transposed.md
%
% Notes:
%   - The old terminal printout is kept.
%   - The Markdown table is much narrower and easier to paste into Word.
%   - Result folders are created automatically.

clc; clear; close all;

%% ===================== 0. Project root and path =====================

ProjectRoot = fileparts(mfilename('fullpath'));
cd(ProjectRoot);
addpath(genpath(ProjectRoot));
rehash;

ResultRoot = fullfile(ProjectRoot, 'Result');
DatasetResultRoot = fullfile(ResultRoot, 'Dataset');

if ~exist(ResultRoot, 'dir')
    mkdir(ResultRoot);
end

if ~exist(DatasetResultRoot, 'dir')
    mkdir(DatasetResultRoot);
end

%% ===================== 1. Diary log =====================

log_file = fullfile(ProjectRoot, 'Overnight_Log.txt');

if exist(log_file, 'file')
    delete(log_file);
end

diary(log_file);
diary on;
cleanupObj = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('\n======================================================\n');
fprintf('  任务开始时间: %s\n', datestr(now));
fprintf('  ProjectRoot: %s\n', ProjectRoot);
fprintf('  ResultRoot : %s\n', DatasetResultRoot);
fprintf('======================================================\n');

%% ===================== 2. 全局配置区域 =====================

datasets    = {'KIN40K','POL','PUMADYN32NM','SARCOS'};
% datasets  = {'SARCOS'};

train_ratio = 0.4;
n_mc        = 10;

% If all are false, this script only reads existing result files and prints/writes summary.
run_log = false;
run_ip  = false;
run_tp  = true;
run_cen = false;
run_nbr = false;

tr_tag = round(train_ratio * 100);

% Pre-create dataset result folders.
for d = 1:numel(datasets)
    SaveFolder = fullfile(DatasetResultRoot, datasets{d});
    if ~exist(SaveFolder, 'dir')
        mkdir(SaveFolder);
    end
end

%% ===================== 3. 方法字典 =====================

methods_dict = {
    'LoG-MOE',      'log_moe_mc%d.mat';
    'LoG-GPOE',     'log_gpoe_mc%d.mat';

    'IP-DAC-MOE',   'moe_tr%d_mc%d.mat';
    'IP-DAC-GPOE',  'gpoe_tr%d_mc%d.mat';
    'IP-DAC-POE',   'poe_tr%d_mc%d.mat';
    'IP-DAC-BCM',   'bcm_tr%d_mc%d.mat';
    'IP-DAC-RBCM',  'rbcm_tr%d_mc%d.mat';

    'IP-AC-MOE',    'moe_ac_tr%d_mc%d.mat';
    'IP-AC-GPOE',   'gpoe_ac_tr%d_mc%d.mat';
    'IP-AC-POE',    'poe_ac_tr%d_mc%d.mat';
    'IP-AC-BCM',    'bcm_ac_tr%d_mc%d.mat';
    'IP-AC-RBCM',   'rbcm_ac_tr%d_mc%d.mat';

    'TP-DAC-MOE',   'moe_tp_tr%d_mc%d.mat';
    'TP-DAC-GPOE',  'gpoe_tp_tr%d_mc%d.mat';
    'TP-DAC-POE',   'poe_tp_tr%d_mc%d.mat';
    'TP-DAC-BCM',   'bcm_tp_tr%d_mc%d.mat';
    'TP-DAC-RBCM',  'rbcm_tp_tr%d_mc%d.mat';

    'TP-AC-MOE',    'moe_ac_tp_tr%d_mc%d.mat';
    'TP-AC-GPOE',   'gpoe_ac_tp_tr%d_mc%d.mat';
    'TP-AC-POE',    'poe_ac_tp_tr%d_mc%d.mat';
    'TP-AC-BCM',    'bcm_ac_tp_tr%d_mc%d.mat';
    'TP-AC-RBCM',   'rbcm_ac_tp_tr%d_mc%d.mat';

    'CEN-MOE',      'moe_cen_tr%d_mc%d.mat';
    'CEN-GPOE',     'gpoe_cen_tr%d_mc%d.mat';
    'CEN-POE',      'poe_cen_tr%d_mc%d.mat';
    'CEN-BCM',      'bcm_cen_tr%d_mc%d.mat';
    'CEN-RBCM',     'rbcm_cen_tr%d_mc%d.mat';

    'NBR-MOE',      'moe_nbr_tr%d_mc%d.mat';
    'NBR-GPOE',     'gpoe_nbr_tr%d_mc%d.mat';
    'NBR-POE',      'poe_nbr_tr%d_mc%d.mat';
    'NBR-BCM',      'bcm_nbr_tr%d_mc%d.mat';
    'NBR-RBCM',     'rbcm_nbr_tr%d_mc%d.mat';
};

methods_names = methods_dict(:,1);
methods_files = methods_dict(:,2);
num_methods   = numel(methods_names);

% Metrics:
%   1 SMSE
%   2 RMSE
%   3 Train time ms/pt
%   4 Test time ms/pt
%   5 Communication train
%   6 Communication test
%   7 Iteration convergence
num_metrics = 7;

mean_results = NaN(numel(datasets), num_methods, num_metrics);
std_results  = NaN(numel(datasets), num_methods, num_metrics);

%% ===================== 4. 运行 & 统计 =====================

for d = 1:numel(datasets)
    dname = datasets{d};
    fprintf('\n>>> 处理数据集: %s <<<\n', dname);

    SaveFolder = fullfile(DatasetResultRoot, dname);
    if ~exist(SaveFolder, 'dir')
        mkdir(SaveFolder);
    end

    %% --------------------- 4.1 Run methods if enabled ---------------------

    for seed = 1:n_mc
        fprintf('\nDataset = %s, seed = %d / %d\n', dname, seed, n_mc);

        if run_log
            try
                ensure_dataset_folder(DatasetResultRoot, dname);
                run_LoGGP_comparison(dname, train_ratio, seed);
            catch ME
                fprintf('\n[ERROR] LoG-GP failed: dataset=%s, seed=%d\n', dname, seed);
                fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
            end
        end

        if run_ip
            try
                ensure_dataset_folder(DatasetResultRoot, dname);
                run_inducingpoint_dataset(dname, 'all', train_ratio, seed);
            catch ME
                fprintf('\n[ERROR] IP methods failed: dataset=%s, seed=%d\n', dname, seed);
                fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
            end
        end

        if run_tp
            try
                ensure_dataset_folder(DatasetResultRoot, dname);
                % num_eval matches NumInducingPoints used in IP-DAC/IP-AC
                % (hardcoded values from run_inducingpoint_dataset.m lines 114-121)
                switch upper(dname)
                    case 'KIN40K',      tp_M = 2100;
                    case 'POL',         tp_M = 2000;
                    case 'PUMADYN32NM', tp_M = 200;
                    case 'SARCOS',      tp_M = 1500;
                    otherwise,          tp_M = 2500;
                end
                fprintf('[TP] Dataset=%s, tp_M=%d\n', dname, tp_M);
                run_testpoint_dataset(dname, 'all', train_ratio, seed, tp_M);
            catch ME
                fprintf('\n[ERROR] TP methods failed: dataset=%s, seed=%d\n', dname, seed);
                fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
            end
        end

        if run_cen
            try
                ensure_dataset_folder(DatasetResultRoot, dname);
                run_centralized_dataset(dname, 'all', train_ratio, seed);
            catch ME
                fprintf('\n[ERROR] Centralized methods failed: dataset=%s, seed=%d\n', dname, seed);
                fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
            end
        end

        if run_nbr
            try
                ensure_dataset_folder(DatasetResultRoot, dname);
                run_neighbor_dataset(dname, 'all', train_ratio, seed);
            catch ME
                fprintf('\n[ERROR] Neighbor methods failed: dataset=%s, seed=%d\n', dname, seed);
                fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
            end
        end
    end

    %% --------------------- 4.2 Collect statistics ---------------------

    fprintf('\nCollecting results for dataset: %s\n', dname);

    for mi = 1:num_methods
        sm   = NaN(1,n_mc);
        rm   = NaN(1,n_mc);
        t_tr = NaN(1,n_mc);
        t_te = NaN(1,n_mc);
        c_tr = NaN(1,n_mc);
        c_te = NaN(1,n_mc);
        it   = NaN(1,n_mc);

        for mc = 1:n_mc
            if contains(methods_files{mi}, 'tr%d')
                fname = sprintf(methods_files{mi}, tr_tag, mc);
            else
                fname = sprintf(methods_files{mi}, mc);
            end

            file_path = fullfile(DatasetResultRoot, dname, fname);

            if exist(file_path, 'file')
                try
                    res = load(file_path);

                    if isfield(res,'smse'), sm(mc) = res.smse; end
                    if isfield(res,'rmse'), rm(mc) = res.rmse; end

                    if isfield(res,'t_train_per_point')
                        t_tr(mc) = res.t_train_per_point;
                    end

                    if isfield(res,'t_test_per_point')
                        t_te(mc) = res.t_test_per_point;
                    end

                    if isfield(res,'comm_train')
                        c_tr(mc) = res.comm_train;
                    else
                        c_tr(mc) = 0;
                    end

                    if isfield(res,'comm_test')
                        c_te(mc) = res.comm_test;
                    else
                        c_te(mc) = 0;
                    end

                    if isfield(res,'iter_converge')
                        it(mc) = res.iter_converge;
                    else
                        it(mc) = NaN;
                    end

                catch ME
                    fprintf('  [读取失败] %s\n', file_path);
                    fprintf('  Reason: %s\n', ME.message);
                end
            else
                fprintf('  [缺失结果] %s\n', file_path);
            end
        end

        mean_results(d,mi,1) = mean(sm,   'omitnan');
        std_results(d,mi,1)  = std(sm,    'omitnan');

        mean_results(d,mi,2) = mean(rm,   'omitnan');
        std_results(d,mi,2)  = std(rm,    'omitnan');

        mean_results(d,mi,3) = mean(t_tr, 'omitnan');
        std_results(d,mi,3)  = std(t_tr,  'omitnan');

        mean_results(d,mi,4) = mean(t_te, 'omitnan');
        std_results(d,mi,4)  = std(t_te,  'omitnan');

        mean_results(d,mi,5) = mean(c_tr, 'omitnan');
        std_results(d,mi,5)  = 0;

        mean_results(d,mi,6) = mean(c_te, 'omitnan');
        std_results(d,mi,6)  = 0;

        mean_results(d,mi,7) = mean(it,   'omitnan');
        std_results(d,mi,7)  = 0;
    end
end

%% ===================== 5. 打印原始宽表到 terminal =====================

metrics_names = {
    'SMSE', ...
    'RMSE', ...
    'Train_T(ms/pt)', ...
    'Test_T(ms/pt)', ...
    'Comm_Train', ...
    'Comm_Test', ...
    'Iter_Converge'
};

agg_list     = {'MOE','GPOE','POE','BCM','RBCM'};
col_prefixes = {'LoG','CEN','IP-DAC','IP-AC','TP-DAC','TP-AC','NBR'};

sep_wide = repmat('=',1,130);
sep_thin = repmat('-',1,130);

for met = 1:num_metrics
    fprintf('\n%s\n', sep_wide);
    fprintf('  Metric: %s  (Train=%.0f%%  MC=%d)\n', metrics_names{met}, train_ratio*100, n_mc);
    fprintf('%s\n', sep_wide);
    fprintf('  Agg    Dataset                   LoG             CEN          IP-DAC           IP-AC          TP-DAC           TP-AC             NBR\n');
    fprintf('  %s\n', sep_thin);

    for a_idx = 1:numel(agg_list)
        agg = agg_list{a_idx};

        for d = 1:numel(datasets)
            if d == 1
                fprintf('  %-4s   %-12s', agg, datasets{d});
            else
                fprintf('         %-12s', datasets{d});
            end

            for col = 1:numel(col_prefixes)
                prefix      = col_prefixes{col};
                target_name = sprintf('%s-%s', prefix, agg);
                mi = find(strcmp(methods_names, target_name), 1);

                if isempty(mi)
                    fprintf('  %14s', '-');
                else
                    mv = mean_results(d,mi,met);
                    sv = std_results(d,mi,met);

                    value_text = format_result_value(mv, sv, met);
                    fprintf('  %14s', value_text);
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

%% ===================== 6. Save summary MAT =====================

SummaryPath = fullfile(DatasetResultRoot, 'All_Methods_Summary.mat');

save(SummaryPath, ...
    'mean_results', 'std_results', ...
    'datasets', 'methods_names', 'metrics_names', ...
    'train_ratio', 'n_mc', 'tr_tag');

fprintf('\nSummary MAT saved to:\n%s\n', SummaryPath);

%% ===================== 7. Save transposed Markdown report =====================

MarkdownPath = fullfile(DatasetResultRoot, 'All_Methods_Summary_transposed.md');

write_transposed_markdown_report( ...
    MarkdownPath, ...
    mean_results, std_results, ...
    datasets, methods_names, ...
    metrics_names, agg_list, col_prefixes, ...
    train_ratio, n_mc);

fprintf('\nTransposed Markdown report saved to:\n%s\n', MarkdownPath);

diary off;

%% ========================================================================
%% Local helper functions
%% ========================================================================

function ensure_dataset_folder(DatasetResultRoot, DatasetName)
    SaveFolder = fullfile(DatasetResultRoot, DatasetName);
    if ~exist(SaveFolder, 'dir')
        mkdir(SaveFolder);
    end
end

function value_text = format_result_value(mv, sv, met)
    if isnan(mv)
        value_text = '-';
        return;
    end

    if met >= 5
        % Communication and iteration columns: integer display.
        % 0 is shown as '-'.
        if mv == 0 || isnan(mv)
            value_text = '-';
        else
            value_text = sprintf('%.0f', mv);
        end
    elseif met == 3 || met == 4
        value_text = sprintf('%.2f±%.2f', mv, sv);
    else
        value_text = sprintf('%.4f±%.4f', mv, sv);
    end
end

function write_transposed_markdown_report( ...
    MarkdownPath, ...
    mean_results, std_results, ...
    datasets, methods_names, ...
    metrics_names, agg_list, col_prefixes, ...
    train_ratio, n_mc)

    fid = fopen(MarkdownPath, 'w');
    if fid < 0
        error('Cannot open Markdown file for writing:\n%s', MarkdownPath);
    end

    fprintf(fid, '# All Methods Summary\n\n');
    fprintf(fid, '**Train ratio:** %.0f%%  \n', train_ratio * 100);
    fprintf(fid, '**Monte Carlo runs:** %d  \n\n', n_mc);
    fprintf(fid, 'Tables are transposed for Word layout. Rows are method families and columns are datasets.\n\n');

    for met = 1:numel(metrics_names)
        metric_name = metrics_names{met};

        % Big metric title
        fprintf(fid, '\n\n# Metric: %s  (Train=%.0f%%, MC=%d)\n\n', ...
            metric_name, train_ratio * 100, n_mc);

        for a_idx = 1:numel(agg_list)
            agg = agg_list{a_idx};

            % Table title
            fprintf(fid, '\n\n## %s - %s\n\n', metric_name, agg);

            % Optional explanatory line
            fprintf(fid, '**Rows:** method families.  **Columns:** datasets.  \n\n');

            % Header
            fprintf(fid, '| Method family |');
            for d = 1:numel(datasets)
                fprintf(fid, ' %s |', datasets{d});
            end
            fprintf(fid, '\n');

            % Alignment row
            fprintf(fid, '|---|');
            for d = 1:numel(datasets)
                fprintf(fid, '---:|');
            end
            fprintf(fid, '\n');

            % Rows: LoG, CEN, IP-DAC, IP-AC, TP-DAC, TP-AC, NBR
            for col = 1:numel(col_prefixes)
                prefix = col_prefixes{col};
                target_name = sprintf('%s-%s', prefix, agg);
                mi = find(strcmp(methods_names, target_name), 1);

                fprintf(fid, '| %s |', prefix);

                for d = 1:numel(datasets)
                    if isempty(mi)
                        value_text = '-';
                    else
                        mv = mean_results(d, mi, met);
                        sv = std_results(d, mi, met);
                        value_text = format_result_value(mv, sv, met);
                    end

                    fprintf(fid, ' %s |', value_text);
                end

                fprintf(fid, '\n');
            end
        end
    end

    fclose(fid);
end