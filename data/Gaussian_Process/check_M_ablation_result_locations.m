%% check_M_ablation_result_locations.m
% Check where the M-ablation result files are stored and whether they contain
% train-time fields needed by plot_inducingpoint_M_tradeoff_pretty_with_time.m.

clear; clc;

Dataset_list = {'KIN40K','POL','PUMADYN32NM','SARCOS'};
Method_list  = {'poe','gpoe','moe','bcm','rbcm'};
M_check      = [500 1000 1500 1800 1900 2000 2100 2500];
seed_check   = 1:3;
train_ratio  = 0.4;
tr_tag       = round(train_ratio * 100);

mainPath = which('run_inducingpoint_dataset_trade_off');
plotPath = which('plot_inducingpoint_M_tradeoff_pretty_with_time');

fprintf('Current pwd:\n%s\n\n', pwd);

if isempty(mainPath)
    warning('run_inducingpoint_dataset_trade_off was not found on MATLAB path.');
    root_main = pwd;
else
    fprintf('run_inducingpoint_dataset_trade_off:\n%s\n\n', mainPath);
    root_main = fileparts(mainPath);
end

if isempty(plotPath)
    fprintf('plot_inducingpoint_M_tradeoff_pretty_with_time not found on path.\n\n');
else
    fprintf('plot_inducingpoint_M_tradeoff_pretty_with_time:\n%s\n\n', plotPath);
end

CandidateRoots = unique({ ...
    pwd, ...
    root_main ...
}, 'stable');

for rr = 1:numel(CandidateRoots)
    root = CandidateRoots{rr};
    fprintf('\n============================================================\n');
    fprintf('Checking root:\n%s\n', root);
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

            for Mi = 1:numel(M_check)
                M = M_check(Mi);

                for si = 1:numel(seed_check)
                    seed = seed_check(si);
                    fname = sprintf('%s_M%d_tr%d_mc%d.mat', method, M, tr_tag, seed);
                    fpath = fullfile(folder, fname);

                    if isfile(fpath)
                        count_found = count_found + 1;
                        S = load(fpath);
                        if isfield(S, 't_train_per_point') || ...
                           (isfield(S, 't_train_total') && isfield(S, 'N_train'))
                            count_with_time = count_with_time + 1;
                        end
                    end
                end
            end

            fprintf('  %-5s files=%3d, with train-time=%3d\n', ...
                upper(method), count_found, count_with_time);
        end
    end
end

fprintf('\nCheck finished.\n');
