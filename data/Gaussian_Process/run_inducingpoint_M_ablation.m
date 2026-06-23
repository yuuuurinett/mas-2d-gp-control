clear; clc;

% Inducing-point ablation: M vs Train time, SMSE, MSLL
DatasetName = 'KIN40K';
train_ratio = 0.4;

% Coarse scan: M from 100 to 2500, step 100 (25 points), across all 5
% aggregation methods. Goal: locate the M-region with the best SMSE/MSLL
% trade-off per method, then (in a separate follow-up run) rescan that
% region with a finer step.
Method_list = {'poe', 'gpoe', 'moe', 'bcm', 'rbcm'};
M_list = 100:100:2500;
seeds = 1:3;

total_runs = numel(Method_list) * numel(M_list) * numel(seeds);
run_counter = 0;

for ci = 1:numel(Method_list)
    CurrentMode = Method_list{ci};

    for mi = 1:numel(M_list)
        M = M_list(mi);
        for si = 1:numel(seeds)
            seed = seeds(si);
            run_counter = run_counter + 1;

            fprintf('\n=================================================\n');
            fprintf('[%d/%d] Dataset: %s | Method: %s | M = %d | seed = %d\n', ...
                run_counter, total_runs, DatasetName, CurrentMode, M, seed);
            fprintf('=================================================\n');

            run_inducingpoint_dataset_trade_off(DatasetName, CurrentMode, train_ratio, seed, M);
        end
    end
end

fprintf('\nAll M-ablation experiments finished.\n');