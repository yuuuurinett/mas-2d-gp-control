%% check_M_ablation_result_locations.m
% Full check for M-ablation result files.
% Expected per dataset/method:
% M = 100:100:2500 -> 25 values
% seeds = 1:3       -> 3 seeds
% total = 75 files per dataset/method

clear; clc;

Dataset_list = {'KIN40K','POL','PUMADYN32NM','SARCOS'};
Method_list  = {'poe','gpoe','moe','bcm','rbcm'};
M_check      = 100:100:2500;
seed_check   = 1:3;
train_ratio  = 0.4;
tr_tag       = round(train_ratio * 100);

mainPath = which('run_inducingpoint_dataset_trade_off');
plotPath = which('plot_inducingpoint_M_tradeoff');

fprintf('Current pwd:\n%s\n\n', pwd);

if isempty(mainPath)
    warning('run_inducingpoint_dataset_trade_off was not found on MATLAB path.');
    root_main = pwd;
else
    fprintf('run_inducingpoint_dataset_trade_off:\n%s\n\n', mainPath);
    root_main = fileparts(mainPath);
end

if isempty(plotPath)
    fprintf('plot_inducingpoint_M_tradeoff not found on path.\n\n');
else
    fprintf('plot_inducingpoint_M_tradeoff:\n%s\n\n', plotPath);
end

CandidateRoots = unique({pwd, root_main}, 'stable');

expected_count = numel(M_check) * numel(seed_check);

for rr = 1:numel(CandidateRoots)
    root = CandidateRoots{rr};

    fprintf('\n============================================================\n');
    fprintf('Checking root:\n%s\n', root);
    fprintf('Expected per dataset/method: %d files\n', expected_count);
    fprintf('============================================================\n');

    for dd = 1:numel(Dataset_list)
        dataset = Dataset_list{dd};
        folder = fullfile(root, 'Result', 'Dataset', dataset);

        fprintf('\nDataset folder: %s\n', folder);

        if ~isfolder(folder)
            fprintf('  [Missing folder]\n');
            continue;
        end

        for mm = 1:numel(Method_list)
            method = Method_list{mm};

            count_found = 0;
            count_with_time = 0;
            count_with_msll = 0;
            count_bad_msll = 0;
            missing_list = {};
            bad_list = {};

            for Mi = 1:numel(M_check)
                M = M_check(Mi);

                for si = 1:numel(seed_check)
                    seed = seed_check(si);

                    fname = sprintf('%s_M%d_tr%d_mc%d.mat', ...
                        method, M, tr_tag, seed);
                    fpath = fullfile(folder, fname);

                    if isfile(fpath)
                        count_found = count_found + 1;

                        S = load(fpath);

                        if isfield(S, 't_train_per_point') || ...
                           (isfield(S, 't_train_total') && isfield(S, 'N_train'))
                            count_with_time = count_with_time + 1;
                        end

                        if isfield(S, 'msll')
                            count_with_msll = count_with_msll + 1;

                            if ~isfinite(S.msll) || abs(S.msll) > 50
                                count_bad_msll = count_bad_msll + 1;
                                bad_list{end+1} = sprintf('%s | MSLL=%.4g', fname, S.msll); %#ok<SAGROW>
                            end
                        else
                            count_bad_msll = count_bad_msll + 1;
                            bad_list{end+1} = sprintf('%s | missing msll', fname); %#ok<SAGROW>
                        end

                    else
                        missing_list{end+1} = fname; %#ok<SAGROW>
                    end
                end
            end

            fprintf('  %-5s files=%3d/%3d, train-time=%3d/%3d, msll=%3d/%3d, bad-msll=%3d\n', ...
                upper(method), ...
                count_found, expected_count, ...
                count_with_time, expected_count, ...
                count_with_msll, expected_count, ...
                count_bad_msll);

            if ~isempty(missing_list)
                fprintf('        Missing examples:\n');
                for k = 1:min(5, numel(missing_list))
                    fprintf('          %s\n', missing_list{k});
                end
                if numel(missing_list) > 5
                    fprintf('          ... %d more missing\n', numel(missing_list)-5);
                end
            end

            if ~isempty(bad_list)
                fprintf('        Bad MSLL examples:\n');
                for k = 1:min(5, numel(bad_list))
                    fprintf('          %s\n', bad_list{k});
                end
                if numel(bad_list) > 5
                    fprintf('          ... %d more bad\n', numel(bad_list)-5);
                end
            end
        end
    end
end

fprintf('\nCheck finished.\n');