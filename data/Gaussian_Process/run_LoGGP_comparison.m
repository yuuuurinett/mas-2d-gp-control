function run_LoGGP_comparison(DatasetName)
%  用 LoG_GP_MultiOutput 在数据集上对比 MOE 和 GPOE 两种聚合方式
%  DatasetName: 'KIN40K' | 'POL' | 'PUMADYN32NM' | 'SARCOS'

rng(0);

%% 1. 加载数据集
fprintf('加载数据集: %s\n', DatasetName);
switch upper(DatasetName)
    case 'KIN40K'
        train = load('KIN40K_train.mat');
        test  = load('KIN40K_test.mat');
        hp    = load('KIN40K_Hyperparameter.mat');
        X_train = train.x;      Y_train = train.y;
        X_test  = test.xtest;   Y_test  = test.ytest;
    case 'POL'
        train = load(fullfile('POL', 'POL_train.mat'));
        test  = load(fullfile('POL', 'POL_test.mat'));
        hp    = load(fullfile('POL', 'POL_Hyperparameter.mat'));
        X_train = train.x;      Y_train = train.y;
        X_test  = test.xtest;   Y_test  = test.ytest;
    case 'PUMADYN32NM'
        train = load(fullfile('PUMADYN32NM', 'PUMADYN32NM_train.mat'));
        test  = load(fullfile('PUMADYN32NM', 'PUMADYN32NM_test.mat'));
        hp    = load(fullfile('PUMADYN32NM', 'PUMADYN32NM_Hyperparameter.mat'));
        X_train = train.x;      Y_train = train.y;
        X_test  = test.xtest;   Y_test  = test.ytest;
    case 'SARCOS'
        train = load(fullfile('SARCOS', 'SARCOS_train.mat'));
        test  = load(fullfile('SARCOS', 'SARCOS_test.mat'));
        hp    = load(fullfile('SARCOS', 'SARCOS_GP_Hyperparameter.mat'));
        % 先用第1组
        hp.SigmaF = hp.SigmaF_set{1};
        hp.SigmaL = hp.SigmaL_set{1};
        hp.SigmaN = hp.SigmaN_set{1};
        X_train = train.sarcos_inv(:, 1:21);
        Y_train = train.sarcos_inv(:, 22:28);
        X_test  = test.sarcos_inv_test(:, 1:21);
        Y_test  = test.sarcos_inv_test(:, 22:28);
    otherwise
        error('未知数据集: %s', DatasetName);
end

[N_train, x_dim] = size(X_train);
y_dim = size(Y_train, 2);
fprintf('Train: %d x %d, Test: %d x %d\n', N_train, x_dim, size(X_test,1), x_dim);
fprintf('SigmaF=%.4f, SigmaN=%.4f\n', hp.SigmaF, hp.SigmaN);

%% 2. 参数设置
Max_LocalGP_DataQuantity = 100;
Max_LocalGP_Quantity     = 100;
methods = {'MOE', 'GPOE'};

%% 3. 训练并评估
N_eval = min(30000, size(X_test,1));
X_eval = X_test(1:N_eval, :);
Y_eval = Y_test(1:N_eval, :);
Y_var_baseline = var(Y_eval);

results = zeros(numel(methods), 5);  % [SMSE, RMSE, NLPD, TrainTime, TestTime]

for method_index = 1:numel(methods)
    method = methods{method_index};
    fprintf('\n--- %s ---\n', method);

    % 初始化 LoG_GP
    log_gp = LoG_GP_MultiOutput( ...
        Max_LocalGP_DataQuantity, Max_LocalGP_Quantity, ...
        x_dim, y_dim, hp.SigmaN, hp.SigmaF, hp.SigmaL);
    log_gp.AggregationMethod = method;

    % 训练
    fprintf('Training...\n');
    tic;
    for n = 1:N_train
        log_gp.update(X_train(n,:)', Y_train(n,:)');
        if mod(n, 1000) == 0
            fprintf('  %d/%d\n', n, N_train);
        end
    end
    t_train = toc;
    fprintf('激活的局部 GP 数量: %d | 耗时: %.2f 秒\n', log_gp.ActivatedGPQuantity, t_train);

    % 评估
    fprintf('Evaluating on %d test points...\n', N_eval);
    mu_pred  = zeros(N_eval, y_dim);
    var_pred = zeros(N_eval, y_dim);
    tic;
    for n = 1:N_eval
        [mu_n, var_n] = log_gp.predict(X_eval(n,:)', Y_eval(n,:)');
        mu_pred(n,:)  = mu_n';
        var_pred(n,:) = var_n';
        if mod(n, 2000) == 0
            fprintf('  %d/%d\n', n, N_eval);
        end
    end
    t_test = toc;
    fprintf('预测耗时: %.2f 秒\n', t_test);

    % 计算指标
    error = Y_eval - mu_pred;
    smse  = mean(error.^2) / Y_var_baseline;
    rmse  = sqrt(mean(error.^2));
    nlpd  = 0.5 * mean(log(2*pi*var_pred) + error.^2 ./ var_pred);

    results(method_index,:) = [mean(smse), mean(rmse), mean(nlpd), t_train, t_test];
end

%% 4. 打印表格
fprintf('\n%s\n', repmat('=',1,72));
fprintf('  %-8s  %10s  %10s  %10s  %10s  %10s\n', ...
    'Method','SMSE','RMSE','NLPD','Train_T(s)','Test_T(s)');
fprintf('  %-8s  %10s  %10s  %10s  %10s  %10s\n', ...
    '------','----','----','----','----------','---------');
for mi = 1:numel(methods)
    fprintf('  %-8s  %10.4f  %10.4f  %10.4f  %10.2f  %10.2f\n', ...
        methods{mi}, results(mi,1), results(mi,2), results(mi,3), results(mi,4), results(mi,5));
end
fprintf('%s\n\n', repmat('=',1,72));

%% 5. 保存和绘图
SaveFolder = fullfile('Result', 'Dataset', DatasetName);
if ~exist(SaveFolder,'dir'), mkdir(SaveFolder); end
fprintf('完成。\n');
end