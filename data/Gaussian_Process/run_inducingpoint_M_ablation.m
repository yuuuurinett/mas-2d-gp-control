%% run_inducingpoint_M_ablation.m
clear; clc;

cd(fileparts(which('run_inducingpoint_dataset_trade_off')));

Dataset_list = {'KIN40K','POL','PUMADYN32NM','SARCOS'};
Method_list  = {'poe','gpoe','moe','bcm','rbcm'};
M_list       = 100:100:2500;
seeds        = 1:3;
train_ratio  = 0.4;
tr_tag       = round(train_ratio * 100);

force_rerun = true;

total_runs = numel(Dataset_list) * numel(Method_list) * numel(M_list) * numel(seeds);
run_counter = 0;
fail_counter = 0;

log_folder = fullfile('Result','Logs');
if ~exist(log_folder, 'dir')
    mkdir(log_folder);
end

log_file = fullfile(log_folder, ['M_ablation_' datestr(now,'yyyymmdd_HHMMSS') '.log']);
diary(log_file);
diary on;

fprintf('============================================================\n');
fprintf('Full M-ablation started at %s\n', datestr(now));
fprintf('Working folder:\n%s\n', pwd);
fprintf('Total runs: %d\n', total_runs);
fprintf('force_rerun = %d\n', force_rerun);
fprintf('============================================================\n');

for di = 1:numel(Dataset_list)
    DatasetName = Dataset_list{di};

    for ci = 1:numel(Method_list)
        CurrentMode = Method_list{ci};

        for mi = 1:numel(M_list)
            M = M_list(mi);

            for si = 1:numel(seeds)
                seed = seeds(si);
                run_counter = run_counter + 1;

                ResultFolder = fullfile('Result','Dataset',DatasetName);
                outFile = fullfile(ResultFolder, sprintf('%s_M%d_tr%d_mc%d.mat', ...
                    CurrentMode, M, tr_tag, seed));

                fprintf('\n============================================================\n');
                fprintf('[%d/%d] Dataset=%s | Method=%s | M=%d | seed=%d\n', ...
                    run_counter, total_runs, DatasetName, CurrentMode, M, seed);
                fprintf('Target file:\n%s\n', outFile);
                fprintf('============================================================\n');

                if ~force_rerun && isfile(outFile)
                    fprintf('[Skip existing]\n');
                    continue;
                end

                try
                    run_inducingpoint_dataset_trade_off( ...
                        DatasetName, CurrentMode, train_ratio, seed, M);

                catch ME
                    fail_counter = fail_counter + 1;

                    fprintf('\n[ERROR]\n');
                    fprintf('Dataset=%s | Method=%s | M=%d | seed=%d\n', ...
                        DatasetName, CurrentMode, M, seed);
                    fprintf('%s\n', ME.message);

                    for kk = 1:numel(ME.stack)
                        fprintf('  at %s line %d\n', ME.stack(kk).name, ME.stack(kk).line);
                    end

                    continue;
                end
            end
        end
    end
end

fprintf('\n============================================================\n');
fprintf('Full M-ablation finished at %s\n', datestr(now));
fprintf('Total runs attempted: %d\n', total_runs);
fprintf('Failures: %d\n', fail_counter);
fprintf('Log file:\n%s\n', log_file);
fprintf('============================================================\n');

diary off;