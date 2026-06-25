%% print_M_train_time_table.m
% Print Low / Medium / High inducing point trade-off table.
% The table includes train time, SMSE, and MSLL.
%
% Recommended use:
%   - Show this table together with M-vs-SMSE/MSLL plots.
%   - Do not use the train-time plot as the only explanation, because timing
%     can be noisy and non-monotonic across seeds/methods.

clear; clc;

ProjectRoot = fileparts(which('run_inducingpoint_dataset_trade_off'));
cd(ProjectRoot);

Dataset_list = {'KIN40K','POL','PUMADYN32NM','SARCOS'};
Method_list  = {'poe','gpoe','moe','bcm','rbcm'};
Method_label = {'POE','GPOE','MOE','BCM','RBCM'};

M_select     = [500 1500 2500];
M_level      = {'Low','Medium','High'};

seeds        = 1:3;
train_ratio  = 0.4;
tr_tag       = round(train_ratio * 100);

ResultRoot = fullfile(ProjectRoot, 'Result', 'Dataset');

OutFolder = fullfile(ProjectRoot, 'Result', 'Figures', 'M_ablation_pretty');
if ~exist(OutFolder, 'dir')
    mkdir(OutFolder);
end

Rows = {};

fprintf('\n============================================================\n');
fprintf('Low / Medium / High M trade-off table\n');
fprintf('Mean over seeds = [%s]\n', num2str(seeds));
fprintf('============================================================\n\n');

for di = 1:numel(Dataset_list)
    dataset = Dataset_list{di};
    fprintf('\nDataset: %s\n', dataset);
    fprintf('%-8s %-7s %-8s | %-15s %-12s %-12s\n', ...
        'Method', 'Level', 'M', 'TrainTime(ms/pt)', 'SMSE', 'MSLL');
    fprintf('%s\n', repmat('-', 1, 72));

    for ci = 1:numel(Method_list)
        method = Method_list{ci};
        method_label = Method_label{ci};

        for mi = 1:numel(M_select)
            M = M_select(mi);
            level = M_level{mi};

            train_time_seed = nan(numel(seeds),1);
            smse_seed       = nan(numel(seeds),1);
            msll_seed       = nan(numel(seeds),1);

            for si = 1:numel(seeds)
                seed = seeds(si);

                fname = sprintf('%s_M%d_tr%d_mc%d.mat', ...
                    method, M, tr_tag, seed);

                fpath = fullfile(ResultRoot, dataset, fname);

                if ~isfile(fpath)
                    fprintf('[Missing] %s\n', fpath);
                    continue;
                end

                S = load(fpath);

                if isfield(S, 'smse')
                    smse_seed(si) = S.smse;
                end

                if isfield(S, 'msll')
                    msll_seed(si) = S.msll;
                end

                if isfield(S, 't_train_per_point')
                    train_time_seed(si) = S.t_train_per_point;
                elseif isfield(S, 't_train_total') && isfield(S, 'N_train') && S.N_train > 0
                    train_time_seed(si) = S.t_train_total / S.N_train * 1000;
                end
            end

            train_time_mean = mean(train_time_seed, 'omitnan');
            train_time_std  = std(train_time_seed, 0, 'omitnan');

            smse_mean = mean(smse_seed, 'omitnan');
            smse_std  = std(smse_seed, 0, 'omitnan');

            msll_mean = mean(msll_seed, 'omitnan');
            msll_std  = std(msll_seed, 0, 'omitnan');

            fprintf('%-8s %-7s M=%-5d | %-15.4f %-12.5g %-12.5g\n', ...
                method_label, level, M, train_time_mean, smse_mean, msll_mean);

            Rows(end+1,:) = { ...
                dataset, method_label, level, M, ...
                train_time_mean, train_time_std, ...
                smse_mean, smse_std, ...
                msll_mean, msll_std ...
            }; %#ok<SAGROW>
        end
    end
end

T = cell2table(Rows, 'VariableNames', { ...
    'Dataset','Method','Level','M', ...
    'TrainTimeMean_ms_per_pt','TrainTimeStd_ms_per_pt', ...
    'SMSE_mean','SMSE_std', ...
    'MSLL_mean','MSLL_std'});

csv_path = fullfile(OutFolder, 'Low_Medium_High_M_tradeoff_table.csv');
writetable(T, csv_path);

fprintf('\nCSV saved to:\n%s\n', csv_path);

%% Export Markdown table
md_path = fullfile(OutFolder, 'Low_Medium_High_M_tradeoff_table.md');
fid = fopen(md_path, 'w');

fprintf(fid, '# Low / Medium / High inducing point trade-off table\n\n');
fprintf(fid, 'Mean over seeds: `%s`\n\n', num2str(seeds));
fprintf(fid, '| Dataset | Method | Level | M | Train Time (ms/pt) | SMSE | MSLL |\n');
fprintf(fid, '|---|---|---:|---:|---:|---:|---:|\n');

 for i = 1:height(T)
    fprintf(fid, '| %s | %s | %s | %d | %.4f | %.5g | %.5g |\n', ...
        T.Dataset{i}, T.Method{i}, T.Level{i}, T.M(i), ...
        T.TrainTimeMean_ms_per_pt(i), ...
        T.SMSE_mean(i), T.MSLL_mean(i));
 end

fclose(fid);

fprintf('Markdown saved to:\n%s\n', md_path);

%% Also print a compact dataset-level average table
fprintf('\n============================================================\n');
fprintf('Compact dataset-level average over methods\n');
fprintf('============================================================\n');

CompactRows = {};

fprintf('%-12s %-7s %-8s | %-15s %-12s %-12s\n', ...
    'Dataset', 'Level', 'M', 'TrainTime(ms/pt)', 'SMSE', 'MSLL');
fprintf('%s\n', repmat('-', 1, 72));

for di = 1:numel(Dataset_list)
    dataset = Dataset_list{di};

    for mi = 1:numel(M_select)
        M = M_select(mi);
        level = M_level{mi};

        mask = strcmp(T.Dataset, dataset) & T.M == M;

        train_time_mean = mean(T.TrainTimeMean_ms_per_pt(mask), 'omitnan');
        smse_mean       = mean(T.SMSE_mean(mask), 'omitnan');
        msll_mean       = mean(T.MSLL_mean(mask), 'omitnan');

        fprintf('%-12s %-7s M=%-5d | %-15.4f %-12.5g %-12.5g\n', ...
            dataset, level, M, train_time_mean, smse_mean, msll_mean);

        CompactRows(end+1,:) = {dataset, level, M, train_time_mean, smse_mean, msll_mean}; %#ok<SAGROW>
    end
end

Tcompact = cell2table(CompactRows, 'VariableNames', { ...
    'Dataset','Level','M','TrainTimeMean_ms_per_pt','SMSE_mean','MSLL_mean'});

compact_csv_path = fullfile(OutFolder, 'Low_Medium_High_M_tradeoff_table_compact.csv');
writetable(Tcompact, compact_csv_path);

fprintf('\nCompact CSV saved to:\n%s\n', compact_csv_path);