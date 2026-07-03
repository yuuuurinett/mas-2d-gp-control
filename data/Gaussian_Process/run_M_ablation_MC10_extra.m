%% run_M_ablation_mixed_MC10.m
% Mixed M-ablation with 10 Monte Carlo runs.
%
% Non-PUMA datasets:
%   M = 100:100:2500
%
% PUMADYN32NM:
%   M = 100:20:500
%
% Results are saved by run_inducingpoint_dataset_trade_off:
%   Result/Dataset/<Dataset>/<method>_M<M>_tr40_mc<seed>.mat

clear; clc; close all;

ProjectRoot = fileparts(mfilename('fullpath'));
cd(ProjectRoot);
addpath(genpath(ProjectRoot));

%% ===================== Configuration =====================

Dataset_list = {'KIN40K','POL','PUMADYN32NM','SARCOS'};
Method_list  = {'poe','gpoe','moe','bcm','rbcm'};

train_ratio = 0.4;
tr_tag = round(train_ratio * 100);

seeds = 1:10;

% Your implementation uses 6 agents in run_inducingpoint_dataset.m
AgentQuantity = 6;

% M settings
M_default = 100:100:2500;
M_puma    = 100:20:500;

ResultRoot = fullfile(ProjectRoot, 'Result', 'Dataset');
if ~exist(ResultRoot, 'dir')
    mkdir(ResultRoot);
end

%% ===================== Diary =====================
diary_file = fullfile(ProjectRoot, 'Overnight_Log_M_ablation_mixed_MC10.txt');

if exist(diary_file, 'file')
    delete(diary_file);
end

diary(diary_file);
diary on;

fprintf('\n============================================================\n');
fprintf('Mixed M-ablation MC10 starts at %s\n', datestr(now));
fprintf('Train ratio = %.0f%%\n', train_ratio * 100);
fprintf('Seeds = %d:%d\n', seeds(1), seeds(end));
fprintf('AgentQuantity = %d\n', AgentQuantity);
fprintf('Default M grid = 100:100:2500\n');
fprintf('PUMA M grid    = 100:20:500\n');
fprintf('============================================================\n');

fprintf('\nAgent/data allocation explanation:\n');
fprintf('  Training data are evenly partitioned across %d agents.\n', AgentQuantity);
fprintf('  Inducing points are globally sampled from the training set.\n');
fprintf('  They are shared as common query locations for all agents.\n');
fprintf('  Error bars will be computed as std over %d Monte Carlo runs.\n', numel(seeds));

%% ===================== Main loop =====================

for di = 1:numel(Dataset_list)
    DatasetName = Dataset_list{di};

    if strcmp(DatasetName, 'PUMADYN32NM')
        M_list = M_puma;
    else
        M_list = M_default;
    end

    SaveFolder = fullfile(ResultRoot, DatasetName);
    if ~exist(SaveFolder, 'dir')
        mkdir(SaveFolder);
    end

    fprintf('\n\n############################################################\n');
    fprintf('Dataset: %s\n', DatasetName);
    fprintf('M list: %s\n', mat2str(M_list));
    fprintf('############################################################\n');

    for mi = 1:numel(M_list)
        M = M_list(mi);

        fprintf('\n============================================================\n');
        fprintf('Dataset=%s | M=%d\n', DatasetName, M);
        fprintf('============================================================\n');

        for si = 1:numel(seeds)
            seed = seeds(si);

            for ai = 1:numel(Method_list)
                method = Method_list{ai};

                save_name = sprintf('%s_M%d_tr%d_mc%d.mat', ...
                    method, M, tr_tag, seed);
                save_path = fullfile(SaveFolder, save_name);

                if isfile(save_path)
                    fprintf('[SKIP] %s\n', save_name);
                    continue;
                end

                fprintf('[RUN] Dataset=%s | Method=%s | M=%d | seed=%d\n', ...
                    DatasetName, upper(method), M, seed);

                try
                    run_inducingpoint_dataset_trade_off( ...
                        DatasetName, method, train_ratio, seed, M);

                catch ME
                    fprintf('[FAILED] Dataset=%s | Method=%s | M=%d | seed=%d\n', ...
                        DatasetName, upper(method), M, seed);
                    fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
                end
            end
        end
    end
end

fprintf('\n============================================================\n');
fprintf('Mixed M-ablation MC10 finished at %s\n', datestr(now));
fprintf('============================================================\n');

diary off;