clear; clc;

% Inducing-point ablation: M vs Train time and SMSE
DatasetName = 'PUMADYN32NM';
CurrentMode = 'poe';
train_ratio = 0.4;

% First quick version. Increase seeds later if needed.
M_list = [500, 1000, 1500, 2000, 2500];
seeds = 1:3;

for mi = 1:numel(M_list)
    M = M_list(mi);
    for si = 1:numel(seeds)
        seed = seeds(si);

        fprintf('\n=================================================\n');
        fprintf('Dataset: %s | Method: %s | M = %d | seed = %d\n', ...
            DatasetName, CurrentMode, M, seed);
        fprintf('=================================================\n');

        run_inducingpoint_dataset_trade_off(DatasetName, CurrentMode, train_ratio, seed, M);
    end
end

fprintf('\nAll M-ablation experiments finished.\n');
