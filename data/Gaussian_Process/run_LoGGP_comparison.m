function run_LoGGP_comparison(DatasetName, train_ratio, seed)
%  用 LoG_GP_MultiOutput 在数据集上对比 MOE 和 GPOE 两种聚合方式
%  DatasetName : 'KIN40K' | 'POL' | 'PUMADYN32NM' | 'SARCOS'
%  train_ratio : 训练集比例，默认 0.4
%  seed        : 随机种子（Monte Carlo 编号），默认 1

if nargin < 2, train_ratio = 0.4; end
if nargin < 3, seed = 1;          end
rng(seed);

%% 1. 加载全量数据
fprintf('加载数据集: %s  (seed=%d, train_ratio=%.0f%%)\n', ...
    DatasetName, seed, train_ratio*100);

switch upper(DatasetName)
    case 'KIN40K'
        tr = load('KIN40K_train.mat');
        te = load('KIN40K_test.mat');
        hp = load('KIN40K_Hyperparameter.mat');
        train_x = tr.x;      train_y = tr.y;
        test_x  = te.xtest;  test_y  = te.ytest;

    case 'POL'
        tr = load(fullfile('POL', 'POL_train.mat'));
        te = load(fullfile('POL', 'POL_test.mat'));
        hp = load(fullfile('POL', 'POL_Hyperparameter.mat'));
        train_x = tr.x;      train_y = tr.y;
        test_x  = te.xtest;  test_y  = te.ytest;

    case 'PUMADYN32NM'
        tr = load(fullfile('PUMADYN32NM', 'PUMADYN32NM_train.mat'));
        te = load(fullfile('PUMADYN32NM', 'PUMADYN32NM_test.mat'));
        hp = load(fullfile('PUMADYN32NM', 'PUMADYN32NM_Hyperparameter.mat'));
        train_x = tr.x;      train_y = tr.y;
        test_x  = te.xtest;  test_y  = te.ytest;

    case 'SARCOS'
        tr = load(fullfile('SARCOS', 'SARCOS_train.mat'));
        te = load(fullfile('SARCOS', 'SARCOS_test.mat'));
        hp = load(fullfile('SARCOS', 'SARCOS_GP_Hyperparameter.mat'));
        hp.SigmaF = hp.SigmaF_set{1};
        hp.SigmaL = hp.SigmaL_set{1};
        hp.SigmaN = hp.SigmaN_set{1};
        train_x = tr.sarcos_inv(:,1:21);       train_y = tr.sarcos_inv(:,22:28);
        test_x  = te.sarcos_inv_test(:,1:21);  test_y  = te.sarcos_inv_test(:,22:28);

    otherwise
        error('未知数据集: %s', DatasetName);
end

%% 2. 与 run_dac_dataset.m 完全相同的数据划分方式
% --- 官方 Train：随机采样 train_ratio% ---
N_official_train = size(train_x, 1);
n_train          = round(N_official_train * train_ratio);
idx_train_all    = randperm(N_official_train);
idx_train_select = idx_train_all(1 : n_train);

X_train = train_x(idx_train_select, :);
Y_train = train_y(idx_train_select, :);

% --- 官方 Test：打乱顺序 ---
N_official_test  = size(test_x, 1);
idx_test_shuffle = randperm(N_official_test);
X_test = test_x(idx_test_shuffle, :);
Y_test = test_y(idx_test_shuffle, :);

[N_train, x_dim] = size(X_train);
y_dim  = size(Y_train, 2);
N_eval = min(30000, size(X_test, 1));
X_eval = X_test(1:N_eval, :);
Y_eval = Y_test(1:N_eval, :);
Y_var_baseline = var(Y_eval);

fprintf('Train(official): %d  Used Train: %d  Test(eval): %d  x_dim: %d  y_dim: %d\n', ...
    N_official_train, N_train, N_eval, x_dim, y_dim);
fprintf('SigmaF=%.4f  SigmaN=%.4f\n', hp.SigmaF, hp.SigmaN);

%% 3. 参数
Max_LocalGP_DataQuantity = 100;
Max_LocalGP_Quantity     = 100;
methods = {'MOE', 'GPOE'};

%% 4. 保存路径
SaveFolder = fullfile('Result', 'Dataset', DatasetName);
if ~exist(SaveFolder, 'dir'), mkdir(SaveFolder); end

results = zeros(numel(methods), 5);  % [SMSE, RMSE, NLPD, Train_T, Test_T]

%% 5. 训练并评估
for method_index = 1:numel(methods)
    method = methods{method_index};
    fprintf('\n--- LoG-%s ---\n', method);

    %% 训练（Train_T 计时开始）
    log_gp = LoG_GP_MultiOutput( ...
        Max_LocalGP_DataQuantity, Max_LocalGP_Quantity, ...
        x_dim, y_dim, hp.SigmaN, hp.SigmaF, hp.SigmaL);
    log_gp.AggregationMethod = method;

    tic;
    for n = 1:N_train
        log_gp.update(X_train(n,:)', Y_train(n,:)');
        if mod(n, 1000) == 0
            fprintf('  训练进度: %d/%d\n', n, N_train);
        end
    end
    t_train = toc;  % Train_T 结束
    fprintf('激活局部GP数量: %d  训练耗时: %.2f 秒\n', log_gp.ActivatedGPQuantity, t_train);

    %% 测试（Test_T 计时开始）
    % LoG-GP 在 test 阶段需要收集各局部GP预测再聚合 → 有通信开销
    mu_pred  = zeros(N_eval, y_dim);
    var_pred = zeros(N_eval, y_dim);

    tic;
    for n = 1:N_eval
        [mu_n, var_n]  = log_gp.predict(X_eval(n,:)', Y_eval(n,:)');
        mu_pred(n, :)  = mu_n';
        var_pred(n, :) = var_n';
        if mod(n, 2000) == 0
            fprintf('  预测进度: %d/%d\n', n, N_eval);
        end
    end
    t_test = toc;  % Test_T 结束
    fprintf('预测耗时: %.2f 秒\n', t_test);

    %% 计算指标
    err  = Y_eval - mu_pred;
    smse = mean(mean(err.^2) ./ Y_var_baseline);
    rmse = mean(sqrt(mean(err.^2)));
    nlpd = mean(mean(0.5*(log(2*pi*var_pred) + err.^2 ./ var_pred)));

    results(method_index, :) = [smse, rmse, nlpd, t_train, t_test];

    fprintf('  SMSE=%.4f  RMSE=%.4f  NLPD=%.4f  Train_T=%.1fs  Test_T=%.1fs\n', ...
        smse, rmse, nlpd, t_train, t_test);

    %% 保存（文件名带 seed 编号）
    method_lower = lower(method);
    save_name    = sprintf('log_%s_mc%d.mat', method_lower, seed);
    save(fullfile(SaveFolder, save_name), ...
        'smse', 'rmse', 'nlpd', 't_train', 't_test', ...
        'method', 'seed', 'train_ratio');
end

%% 6. 打印对比表格
fprintf('\n%s\n', repmat('=', 1, 72));
fprintf('  %-12s  %8s  %8s  %8s  %10s  %10s\n', ...
    'Method', 'SMSE', 'RMSE', 'NLPD', 'Train_T(s)', 'Test_T(s)');
fprintf('  %s\n', repmat('-', 1, 68));
method_labels = {'LoG-MOE', 'LoG-GPOE'};
for mi = 1:numel(methods)
    fprintf('  %-12s  %8.4f  %8.4f  %8.4f  %10.2f  %10.2f\n', ...
        method_labels{mi}, results(mi,1), results(mi,2), ...
        results(mi,3), results(mi,4), results(mi,5));
end
fprintf('%s\n\n', repmat('=', 1, 72));

fprintf('[%s] seed=%d 完成。\n\n', DatasetName, seed);
end