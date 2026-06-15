clear; clc; close all;

DatasetName = 'SARCOS';
Method = 'poe';
train_ratio = 0.4;
tr_tag = round(train_ratio * 100);

M_list = [500, 1000, 1500, 2000, 2500];
seeds = 1:3;

% This matches the SaveFolder in run_inducingpoint_dataset.m:
% SaveFolder = fullfile('Result','Dataset',DatasetName);
ResultFolder = fullfile('Result','Dataset',DatasetName);

SMSE_all   = nan(numel(M_list), numel(seeds));
TrainT_all = nan(numel(M_list), numel(seeds));
TestT_all  = nan(numel(M_list), numel(seeds));
KernelT_all = nan(numel(M_list), numel(seeds));
MeanT_all   = nan(numel(M_list), numel(seeds));

for mi = 1:numel(M_list)
    M = M_list(mi);

    for si = 1:numel(seeds)
        seed = seeds(si);

        file_name = sprintf('%s_M%d_tr%d_mc%d.mat', Method, M, tr_tag, seed);
        file_path = fullfile(ResultFolder, file_name);

        if ~exist(file_path, 'file')
            warning('Missing file: %s', file_path);
            continue;
        end

        S = load(file_path);

        if isfield(S, 'smse')
            SMSE_all(mi, si) = S.smse;
        end

        if isfield(S, 't_train_per_point')
            TrainT_all(mi, si) = S.t_train_per_point;   % ms / training point
        elseif isfield(S, 't_train_total') && isfield(S, 'N_train')
            TrainT_all(mi, si) = S.t_train_total / S.N_train * 1000;
        end

        if isfield(S, 't_test_per_point')
            TestT_all(mi, si) = S.t_test_per_point;     % ms / test point
        elseif isfield(S, 't_test_total') && isfield(S, 'N_eval')
            TestT_all(mi, si) = S.t_test_total / S.N_eval * 1000;
        end

        % Optional detailed test-time components, normalized as ms/test point
        if isfield(S, 't_ip_test_kernel') && isfield(S, 'N_eval')
            KernelT_all(mi, si) = S.t_ip_test_kernel / S.N_eval * 1000;
        end
        if isfield(S, 't_ip_test_mean') && isfield(S, 'N_eval')
            MeanT_all(mi, si) = S.t_ip_test_mean / S.N_eval * 1000;
        end
    end
end

SMSE_mean   = mean(SMSE_all,   2, 'omitnan');
SMSE_std    = std(SMSE_all,  0, 2, 'omitnan');
TrainT_mean = mean(TrainT_all, 2, 'omitnan');
TrainT_std  = std(TrainT_all,0, 2, 'omitnan');
TestT_mean  = mean(TestT_all,  2, 'omitnan');
TestT_std   = std(TestT_all, 0, 2, 'omitnan');

fprintf('\n============================================================\n');
fprintf('Inducing-point ablation: %s / %s\n', DatasetName, upper(Method));
fprintf('============================================================\n');
fprintf('%10s  %20s  %18s  %18s\n', ...
    'M', 'Train_T(ms/pt)', 'Test_T(ms/pt)', 'SMSE');
fprintf('------------------------------------------------------------\n');

for mi = 1:numel(M_list)
    fprintf('%10d  %9.4f ± %-8.4f  %8.4f ± %-8.4f  %8.5f ± %-8.5f\n', ...
        M_list(mi), ...
        TrainT_mean(mi), TrainT_std(mi), ...
        TestT_mean(mi), TestT_std(mi), ...
        SMSE_mean(mi), SMSE_std(mi));
end

%% Figure 1: M vs training time, log y-axis
figure;
errorbar(M_list, TrainT_mean, TrainT_std, '-o', ...
    'LineWidth', 1.5, 'MarkerSize', 7);
set(gca, 'YScale', 'log');
grid on;
xlabel('Number of inducing points M');
ylabel('Training time (ms / training point, log scale)');
title(sprintf('%s / %s: M vs training time', DatasetName, upper(Method)));

%% Figure 2: M vs SMSE
figure;
errorbar(M_list, SMSE_mean, SMSE_std, '-o', ...
    'LineWidth', 1.5, 'MarkerSize', 7);
grid on;
xlabel('Number of inducing points M');
ylabel('SMSE');
title(sprintf('%s / %s: M vs SMSE', DatasetName, upper(Method)));

%% Figure 3: runtime-performance trade-off
figure;
plot(TrainT_mean, SMSE_mean, '-o', ...
    'LineWidth', 1.5, 'MarkerSize', 7);
grid on;
xlabel('Training time (ms / training point)');
ylabel('SMSE');
title(sprintf('%s / %s: runtime-performance trade-off', DatasetName, upper(Method)));

for mi = 1:numel(M_list)
    text(TrainT_mean(mi), SMSE_mean(mi), sprintf('  M=%d', M_list(mi)), ...
        'FontSize', 10);
end

%% Optional Figure 4: M vs test time
figure;
errorbar(M_list, TestT_mean, TestT_std, '-o', ...
    'LineWidth', 1.5, 'MarkerSize', 7);
set(gca, 'YScale', 'log');
grid on;
xlabel('Number of inducing points M');
ylabel('Test time (ms / test point, log scale)');
title(sprintf('%s / %s: M vs test time', DatasetName, upper(Method)));
