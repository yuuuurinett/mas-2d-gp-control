function run_neighbor_dataset(DatasetName, CurrentMode, train_ratio, seed)

if nargin < 3, train_ratio = 0.4; end
if nargin < 4, seed = 1;          end
rng(seed);

fprintf('\n[邻域静态聚合] %s  seed=%d  tr=%.0f%%\n', DatasetName, seed, train_ratio*100);

%% 1. 加载数据集
switch upper(DatasetName)
    case 'KIN40K'
        tr=load('KIN40K_train.mat'); te=load('KIN40K_test.mat');
        hp=load('KIN40K_Hyperparameter.mat');
        train_x=tr.x; train_y=tr.y; test_x=te.xtest; test_y=te.ytest;
    case 'POL'
        tr=load(fullfile('POL','POL_train.mat')); te=load(fullfile('POL','POL_test.mat'));
        hp=load(fullfile('POL','POL_Hyperparameter.mat'));
        train_x=tr.x; train_y=tr.y; test_x=te.xtest; test_y=te.ytest;
    case 'PUMADYN32NM'
        tr=load(fullfile('PUMADYN32NM','PUMADYN32NM_train.mat')); 
        te=load(fullfile('PUMADYN32NM','PUMADYN32NM_test.mat'));
        hp=load(fullfile('PUMADYN32NM','PUMADYN32NM_Hyperparameter.mat'));
        train_x=tr.x; train_y=tr.y; test_x=te.xtest; test_y=te.ytest;
    case 'SARCOS'
        tr = load(fullfile('SARCOS','SARCOS_train.mat')); 
        te = load(fullfile('SARCOS','SARCOS_test.mat'));
        hp_raw = load(fullfile('SARCOS','SARCOS_GP_Hyperparameter.mat'));
        hp.SigmaF = mean(cell2mat(hp_raw.SigmaF_set)); 
        hp.SigmaN = mean(cell2mat(hp_raw.SigmaN_set)); 
        hp.SigmaL = mean(cell2mat(hp_raw.SigmaL_set'), 2);
        train_x = tr.sarcos_inv(:,1:21);       train_y = tr.sarcos_inv(:,22:28); 
        test_x  = te.sarcos_inv_test(:,1:21);  test_y  = te.sarcos_inv_test(:,22:28);
end

%% 2. 超参数提取 & 归一化
if size(hp.SigmaL,1)>1 && size(hp.SigmaL,2)>1, hp.SigmaL = mean(hp.SigmaL,1); end
SigmaL = hp.SigmaL(:);
if numel(hp.SigmaF)>1, hp.SigmaF = mean(hp.SigmaF); end
if numel(hp.SigmaN)>1, hp.SigmaN = mean(hp.SigmaN); end
SigmaF = hp.SigmaF;  SigmaN = hp.SigmaN;

[N_total, x_dim] = size(train_x);
y_dim  = size(train_y, 2);
n_train = round(N_total * train_ratio);
idx_tr  = randperm(N_total, n_train);
X_train = train_x(idx_tr, :);
Y_train = train_y(idx_tr, :);

% 测试集洗牌
idx_te  = randperm(size(test_x, 1));
X_test  = test_x(idx_te, :);
Y_test  = test_y(idx_te, :);

X_mean = mean(X_train, 1);
X_std  = std(X_train, 0, 1);  X_std(X_std==0) = 1;
if ~(max(abs(X_mean))<1e-2 && max(abs(X_std-1))<1e-2)
    X_train = (X_train - X_mean) ./ X_std;
    X_test  = (X_test  - X_mean) ./ X_std;
    SigmaL  = SigmaL(:) ./ X_std(:);
end

Y_mean = mean(Y_train, 1);
Y_std  = std(Y_train, 0, 1);  Y_std(Y_std==0) = 1;
if ~(max(abs(Y_mean))<1e-2 && max(abs(Y_std-1))<1e-2)
    Y_train = (Y_train - Y_mean) ./ Y_std;
    mean_Y_std = mean(Y_std);
    SigmaF  = SigmaF / mean_Y_std;
    SigmaN  = SigmaN / mean_Y_std;
else
    Y_mean = zeros(1, y_dim);
    Y_std  = ones(1, y_dim);
end
prior_var = SigmaF^2;

N_eval     = min(3000, size(X_test, 1));
X_eval     = X_test(1:N_eval, :);
Y_eval_raw = Y_test(1:N_eval, :);       
Y_eval     = (Y_eval_raw - Y_mean) ./ Y_std;  
Y_var_base = var(Y_eval_raw, 0, 1);    

%% 3. 拓扑结构与训练局部专家
AgentQuantity = 6;
MAS = Manipulator_2D_2DoF_SetMASTopology(AgentQuantity, 1);
L = MAS.Agent_Topology.LaplacianMatrix;

MaxData = floor(n_train / AgentQuantity);
LocalGP_set = cell(AgentQuantity, 1);
for n = 1:AgentQuantity
    idx = (n-1)*MaxData+1 : n*MaxData;
    LocalGP_set{n} = LocalGP_MultiOutput(x_dim, y_dim, MaxData, SigmaN, SigmaF, SigmaL);
    LocalGP_set{n}.add_Alldata(X_train(idx, :)', Y_train(idx, :)');
end

%% =========================================================================
%% 4. 统一预计算阶段 (极其耗时，提出到最外层)
%% =========================================================================
if strcmpi(CurrentMode, 'all')
    method_list = {'moe', 'poe', 'gpoe', 'bcm', 'rbcm'};
else
    method_list = {lower(CurrentMode)};
end

fprintf('正在生成 %d 个测试点的局部推断...\n', N_eval);
tic;
mu_all  = zeros(AgentQuantity, y_dim, N_eval);
var_all = zeros(AgentQuantity, y_dim, N_eval);
for n = 1:AgentQuantity
    for t = 1:N_eval
        [m, v] = LocalGP_set{n}.predict(X_eval(t, :)');
        mu_all(n, :, t) = m';
        var_all(n, :, t) = v';
    end
end
Precompute_Time = toc;
fprintf('预计算完成！耗时: %.2f 秒\n', Precompute_Time);

%% =========================================================================
%% 5. 执行静态邻域聚合 (3D 矩阵向量化加速)
%% =========================================================================
for mi = 1:numel(method_list)
    cur_m = method_list{mi};
    fprintf('  处理方法: %s\n', upper(cur_m));
    tic; 
    
    agent_smse = zeros(AgentQuantity, 1);
    agent_rmse = zeros(AgentQuantity, 1);
    agent_nlpd = zeros(AgentQuantity, 1);
    
    for n = 1:AgentQuantity
        % 1. 明确邻域集合
        Neighbor_Indices = [n, find(L(n, :) < 0)];  
        Num_Neighbors = numel(Neighbor_Indices);    
        
        % 2. 批量提取该智能体整个邻域的所有预测 [Num_Neighbors, y_dim, N_eval]
        N_Means = mu_all(Neighbor_Indices, :, :);  
        N_Vars  = var_all(Neighbor_Indices, :, :); 
        
        % 3. 批量计算动态权重 Beta
        Beta = max(0.5 * (log(prior_var) - log(N_Vars)), eps); 
        
        % 4. 核心聚合公式 (消灭 for t 和 for d 循环，沿第 1 维度/邻居维 求和)
        switch cur_m
            case 'moe'
                M_Agg = mean(N_Means, 1);
                V_Agg = mean(N_Vars + N_Means.^2, 1) - M_Agg.^2;
            case 'poe'
                Sum_Prec = sum(1 ./ N_Vars, 1);
                M_Agg = sum(N_Means ./ N_Vars, 1) ./ Sum_Prec;
                V_Agg = 1 ./ Sum_Prec;
            case 'gpoe'
                Sum_Prec = sum(Beta ./ N_Vars, 1);
                M_Agg = sum(Beta .* N_Means ./ N_Vars, 1) ./ Sum_Prec;
                V_Agg = 1 ./ Sum_Prec;
            case 'bcm'
                Sum_Prec = max(sum(1 ./ N_Vars, 1) - (Num_Neighbors - 1) / prior_var, eps);
                M_Agg = sum(N_Means ./ N_Vars, 1) ./ Sum_Prec;
                V_Agg = 1 ./ Sum_Prec;
            case 'rbcm'
                Sum_Prec = max(sum(Beta ./ N_Vars, 1) + (1 - sum(Beta, 1)) / prior_var, eps);
                M_Agg = sum(Beta .* N_Means ./ N_Vars, 1) ./ Sum_Prec;
                V_Agg = 1 ./ Sum_Prec;
        end
        
        % 5. 调整维度从 [1, y_dim, N_eval] 变为 [N_eval, y_dim]
        Mean_Aggregated = permute(M_Agg, [3, 2, 1]);
        Var_Aggregated  = permute(V_Agg, [3, 2, 1]);
        Var_Aggregated  = max(Var_Aggregated, SigmaN^2); % 下界保护
        
        % 6. 反归一化预测均值，并在原始空间计算误差
        Mean_Aggregated_orig = Mean_Aggregated .* repmat(Y_std, N_eval, 1) + repmat(Y_mean, N_eval, 1);
        Var_orig = Var_Aggregated .* repmat(Y_std.^2, N_eval, 1);
        err = Y_eval_raw - Mean_Aggregated_orig;
        
        % 7. 记录该 Agent 视角的评估指标
        agent_smse(n) = mean(mean(err.^2) ./ Y_var_base);
        agent_rmse(n) = mean(sqrt(mean(err.^2)));
        agent_nlpd(n) = mean(mean(0.5*(log(2*pi*Var_orig) + err.^2./Var_orig)));
    end
    
    % 取全网 6 个智能体视角的平均值作为方法最终得分
    smse = mean(agent_smse);
    rmse = mean(agent_rmse);
    nlpd = mean(agent_nlpd);
    
    % [公平时间计算]: 均摊预计算时间 + 本次矩阵融合运算时间
    t_test = Precompute_Time/numel(method_list) + toc; 
    t_train = 0; % NBR 没有离线共识时间

    fprintf('    SMSE=%.4f  RMSE=%.4f  NLPD=%.4f  Time=%.1fs\n', smse, rmse, nlpd, t_test);

    % 保存结果
    SaveFolder = fullfile('Result', 'Dataset', DatasetName);
    if ~exist(SaveFolder,'dir'), mkdir(SaveFolder); end
    save(fullfile(SaveFolder, sprintf('%s_nbr_mc%d.mat', cur_m, seed)), ...
        'smse', 'rmse', 'nlpd', 't_train', 't_test', 'seed', 'train_ratio', 'cur_m');
end
fprintf('\n[%s] 邻域聚合 done. seed=%d tr=%d%%\n\n', DatasetName, seed, round(train_ratio*100));
end