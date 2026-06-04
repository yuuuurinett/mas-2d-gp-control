function run_inducingpoint_dataset(DatasetName, CurrentMode, train_ratio, seed)
%  诱导点聚合方法（IP-DAC / IP-AC）
if nargin < 3, train_ratio = 0.4; end
if nargin < 4, seed = 1;          end
rng(seed);

fprintf('\n[诱导点] %s  seed=%d  tr=%.0f%%\n', DatasetName, seed, train_ratio*100);

%% 1. 加载数据集
switch upper(DatasetName)
    case 'KIN40K'
        tr = load('KIN40K_train.mat');          te = load('KIN40K_test.mat');
        hp = load('KIN40K_Hyperparameter.mat');
        train_x = tr.x;  train_y = tr.y;  test_x = te.xtest;  test_y = te.ytest;
    case 'POL'
        tr = load(fullfile('POL','POL_train.mat'));
        te = load(fullfile('POL','POL_test.mat'));
        hp = load(fullfile('POL','POL_Hyperparameter.mat'));
        train_x = tr.x;  train_y = tr.y;  test_x = te.xtest;  test_y = te.ytest;
    case 'PUMADYN32NM'
        tr = load(fullfile('PUMADYN32NM','PUMADYN32NM_train.mat'));
        te = load(fullfile('PUMADYN32NM','PUMADYN32NM_test.mat'));
        hp = load(fullfile('PUMADYN32NM','PUMADYN32NM_Hyperparameter.mat'));
        train_x = tr.x;  train_y = tr.y;  test_x = te.xtest;  test_y = te.ytest;
    case 'SARCOS'
        tr  = load(fullfile('SARCOS','SARCOS_train.mat'));
        te  = load(fullfile('SARCOS','SARCOS_test.mat'));
        hp_raw = load(fullfile('SARCOS','SARCOS_GP_Hyperparameter.mat'));
        hp.SigmaF = mean(cell2mat(hp_raw.SigmaF_set));
        hp.SigmaN = mean(cell2mat(hp_raw.SigmaN_set));
        hp.SigmaL = mean(cell2mat(hp_raw.SigmaL_set'), 2);
        train_x = tr.sarcos_inv(:,1:21);       train_y = tr.sarcos_inv(:,22:28);
        test_x  = te.sarcos_inv_test(:,1:21);  test_y  = te.sarcos_inv_test(:,22:28);
    otherwise,  error('未知数据集: %s', DatasetName);
end

%% 2. 归一化
if size(hp.SigmaL,1)>1 && size(hp.SigmaL,2)>1, hp.SigmaL = mean(hp.SigmaL,1); end
SigmaL = hp.SigmaL(:);
if numel(hp.SigmaF)>1, hp.SigmaF = mean(hp.SigmaF); end
if numel(hp.SigmaN)>1, hp.SigmaN = mean(hp.SigmaN); end
SigmaF = hp.SigmaF;  SigmaN = hp.SigmaN;

num_train_samples = round(size(train_x,1) * train_ratio);
train_indices     = randperm(size(train_x,1), num_train_samples);
X_train = train_x(train_indices,:);  Y_train = train_y(train_indices,:);

test_indices = randperm(size(test_x,1));
X_test  = test_x(test_indices,:);    Y_test  = test_y(test_indices,:);

X_mean = mean(X_train,1);  X_std = std(X_train,0,1);  X_std(X_std==0) = 1;
if ~(max(abs(X_mean))<1e-2 && max(abs(X_std-1))<1e-2)
    X_train = (X_train - X_mean) ./ X_std;
    X_test  = (X_test  - X_mean) ./ X_std;
    SigmaL  = SigmaL ./ X_std(:);
end
Y_mean = mean(Y_train,1);  Y_std = std(Y_train,0,1);  Y_std(Y_std==0) = 1;
if max(abs(Y_mean))<1e-2 && max(abs(Y_std-1))<1e-2
    Y_mean = zeros(1,size(Y_train,2));  Y_std = ones(1,size(Y_train,2));
else
    Y_train = (Y_train - Y_mean) ./ Y_std;
    SigmaF  = SigmaF / mean(Y_std);
    SigmaN  = SigmaN / mean(Y_std);
end
prior_var = SigmaF^2;

[num_train, input_dim] = size(X_train);
output_dim = size(Y_train,2);
num_eval   = min(3000, size(X_test,1));
X_eval     = X_test(1:num_eval,:);
Y_eval     = Y_test(1:num_eval,:);
Y_var_base = var(Y_eval, 0, 1);
fprintf('Train=%d  Test=%d  x_dim=%d  y_dim=%d\n', num_train, num_eval, input_dim, output_dim);

%% 3. 分布式系统参数
num_agents        = 6;
dac_step_size     = 0.01;  
dac_gain          = 10;     
max_data_per_agent = min(floor(num_train / num_agents), 3000);

switch upper(DatasetName)
    case {'SARCOS','POL'},  num_inducing_points = 2500;
    otherwise,              num_inducing_points = 2000;
end

MultiAgentSystem = Manipulator_2D_2DoF_SetMASTopology(num_agents, 1);

Laplacian = MultiAgentSystem.Agent_Topology.LaplacianMatrix;

% 各 agent 的邻居数量 |N_i|（用于 ET 阈值计算）
neighbor_count_per_agent = sum(Laplacian < 0, 2);  % [num_agents × 1]
max_neighbor_count       = max(neighbor_count_per_agent);


% a_param：满足论文约束 0 < a < 1/|N_i|，取 0.9/max|N_i| 保证所有 agent 满足约束
et_sigma = 0.2;
et_a     = 0.9 / max_neighbor_count;
fprintf('ET 参数: sigma=%.2f  a=%.4f  (max|N_i|=%d，约束上界=%.4f)\n', ...
    et_sigma, et_a, max_neighbor_count, 1/max_neighbor_count);

%% 4. 训练局部 GP
t_start = tic;
local_gp_set = cell(num_agents, 1);
for agent_idx = 1:num_agents
    data_idx = (agent_idx-1)*max_data_per_agent+1 : min(agent_idx*max_data_per_agent, num_train);
    local_gp_set{agent_idx} = LocalGP_MultiOutput(input_dim, output_dim, ...
        max_data_per_agent, SigmaN, SigmaF, SigmaL);
    local_gp_set{agent_idx}.add_Alldata(X_train(data_idx,:)', Y_train(data_idx,:)');
    local_gp_set{agent_idx}.tau   = 1e-8;
    local_gp_set{agent_idx}.delta = 0.01;
end
t_train_gp = toc(t_start);
fprintf('局部GP训练完成: %.4fs\n', t_train_gp);

% 随机选取诱导点坐标（从训练集中抽取）
inducing_indices     = randperm(num_train, num_inducing_points);
InducingPoints_X     = X_train(inducing_indices,:)';  % [input_dim × num_inducing_points]

%% 5. 方法列表
dac_methods = {'poe','gpoe','moe','bcm','rbcm'};
ac_methods  = {'moe_ac','gpoe_ac','poe_ac','bcm_ac','rbcm_ac'};
if strcmpi(CurrentMode,'all')
    AllModes = [dac_methods, ac_methods];
else
    AllModes = {lower(CurrentMode)};
end

%% 6. 批量预计算每个局部 GP 在所有诱导点上的预测
% P_matrix_xxx：信息向量矩阵，维度 [2*output_dim × num_agents × num_inducing_points]
%   奇数行 2d-1：精度加权均值项（numerator）
%   偶数行 2d  ：精度项（denominator）
p_dim = 2 * output_dim;
P_matrix_poe  = zeros(p_dim, num_agents, num_inducing_points);
P_matrix_gpoe = zeros(p_dim, num_agents, num_inducing_points);
P_matrix_moe  = zeros(p_dim, num_agents, num_inducing_points);
P_matrix_bcm  = zeros(p_dim, num_agents, num_inducing_points);
P_matrix_rbcm = zeros(p_dim, num_agents, num_inducing_points);

% 同时保存局部预测均值和方差
local_mu_at_inducing  = zeros(num_agents, output_dim, num_inducing_points);
local_var_at_inducing = zeros(num_agents, output_dim, num_inducing_points);

fprintf('[批量预计算] 正在计算 %d 个诱导点的局部预测...\n', num_inducing_points);
tic;
for agent_idx = 1:num_agents
    % batch_predict_external 一次性对所有诱导点做预测
    % 输出: mu_batch [num_inducing_points × output_dim], var_batch [num_inducing_points × 1]
    [mu_batch, var_batch] = batch_predict_external( ...
        local_gp_set{agent_idx}, InducingPoints_X, SigmaN, SigmaF);

    for inducing_idx = 1:num_inducing_points
        local_mu_i  = mu_batch(inducing_idx,:)';         % [output_dim × 1] 局部预测均值
        local_var_i = repmat(var_batch(inducing_idx), output_dim, 1); % [output_dim × 1] 局部预测方差

        local_mu_at_inducing(agent_idx,:,inducing_idx)  = local_mu_i';
        local_var_at_inducing(agent_idx,:,inducing_idx) = local_var_i';

        for dim_idx = 1:output_dim
            % 对方差加下界，防止除以极小值导致数值爆炸
            safe_var = max(local_var_i(dim_idx), SigmaN^2);

            % gPoE 的权重 beta_i = 0.5*(log(prior_var) - log(var_i))
            % 截断到 [eps, 10] 防止数值问题（RBCM 对 beta 敏感）
            beta_gpoe = max(min(0.5*(log(prior_var) - log(safe_var)), 10), eps);

            % --- PoE 信息向量 ---
            P_matrix_poe(2*dim_idx-1, agent_idx, inducing_idx) = num_agents * local_mu_i(dim_idx) / safe_var;
            P_matrix_poe(2*dim_idx,   agent_idx, inducing_idx) = num_agents / safe_var;

            % --- gPoE 信息向量 ---
            P_matrix_gpoe(2*dim_idx-1, agent_idx, inducing_idx) = num_agents * beta_gpoe * local_mu_i(dim_idx) / safe_var;
            P_matrix_gpoe(2*dim_idx,   agent_idx, inducing_idx) = num_agents * beta_gpoe / safe_var;

            % --- MoE 信息向量 ---
            P_matrix_moe(2*dim_idx-1, agent_idx, inducing_idx) = num_agents * local_mu_i(dim_idx);
            P_matrix_moe(2*dim_idx,   agent_idx, inducing_idx) = num_agents * (safe_var + local_mu_i(dim_idx)^2);

            % --- BCM 信息向量（含先验精度修正项）---
            P_matrix_bcm(2*dim_idx-1, agent_idx, inducing_idx) = num_agents * local_mu_i(dim_idx) / safe_var;
            P_matrix_bcm(2*dim_idx,   agent_idx, inducing_idx) = num_agents / safe_var - (num_agents-1) / prior_var;

            % --- rBCM 信息向量（beta 加权 + 先验修正）---
            P_matrix_rbcm(2*dim_idx-1, agent_idx, inducing_idx) = num_agents * beta_gpoe * local_mu_i(dim_idx) / safe_var;
            P_matrix_rbcm(2*dim_idx,   agent_idx, inducing_idx) = num_agents * beta_gpoe / safe_var + (1 - num_agents*beta_gpoe) / prior_var;
        end
    end
end
fprintf('预计算完成，耗时: %.2fs\n', toc);

SaveFolder = fullfile('Result','Dataset',DatasetName);
if ~exist(SaveFolder,'dir'), mkdir(SaveFolder); end
tr_tag = round(train_ratio * 100);

%% 7. 各方法主循环
for method_idx = 1:numel(AllModes)
    current_method = AllModes{method_idx};
    fprintf('\n[%d/%d] 方法: %s\n', method_idx, numel(AllModes), current_method);
    tic;

    switch lower(current_method)
       case dac_methods
            %% --- IP-DAC：诱导点分布式平均共识（DAC）---
            switch lower(current_method)
                case 'gpoe', Pi = P_matrix_gpoe;
                case 'poe',  Pi = P_matrix_poe;
                case 'bcm',  Pi = P_matrix_bcm;
                case 'rbcm', Pi = P_matrix_rbcm;
                otherwise,   Pi = P_matrix_moe;
            end

            % =======================================================
            % 核心修复 1: 纯正 DAC 初始化 (直接以 Pi 为初始状态)
            % =======================================================
            Xi              = Pi;  % 真实的连续时间状态
            Xi_last_trigger = Pi;  % 广播快照状态 (用于邻居间通信)
            
            % 核心修复 2: 独立通信计数器，大小为 [num_agents x num_inducing_points]
            % 精确记录“每个智能体”对“每个单独诱导点”的触发次数
            trigger_count_per_agent = zeros(num_agents, num_inducing_points);

            dac_iter  = 0;
            max_iters = 3000;

            % 导师建议的 ET 参数与防 Zeno Bound
            c_0 = 0.1;           % Bound 的初始大小 (可微调)
            alpha_decay = 0.02;  % Bound 的指数衰减率 (可微调)

            while dac_iter < max_iters
                dac_iter  = dac_iter + 1;
                Xi_prev   = Xi;  

                %% Step 1：计算邻居通过网络传来的驱动力 (L * Xi_last_trigger)
                % 【注意】这里必须严格使用邻居的“广播状态” Xi_last_trigger，禁止“作弊”使用真实 Xi
                L_Xi = zeros(size(Xi_last_trigger));
                for agent_idx = 1:num_agents
                    % reshape Laplacian 变三维以便于矩阵广播相乘
                    L_Xi(:, agent_idx, :) = sum(Xi_last_trigger .* reshape(Laplacian(agent_idx,:), 1, num_agents, 1), 2);
                end

                %% Step 2：真实状态演化 (遵循标准 DAC 动力学: \dot{X} = -L * \hat{X})
                for agent_idx = 1:num_agents
                    Xi(:, agent_idx, :) = Xi(:, agent_idx, :) - dac_step_size * dac_gain * L_Xi(:, agent_idx, :);
                end

                %% Step 3：ET 触发判定 (Point-wise 完全解耦，各诱导点独立触发)
                % 引入防 Zeno 现象的指数衰减 Bound
                Bound = c_0 * exp(-alpha_decay * dac_iter);
                
                for agent_idx = 1:num_agents
                    % 当前 agent 的误差矩阵 [p_dim x 1 x num_inducing_points]
                    error_i = Xi_last_trigger(:, agent_idx, :) - Xi(:, agent_idx, :);
                    
                    % 【核心修复 3】在特征维度(p_dim)上求平方和，千万不能用 (:)
                    % 结果得到大小为 [num_inducing_points x 1] 的向量，保留诱导点的独立性
                    norm_e_sq = squeeze(sum(error_i.^2, 1));     
                    norm_z_sq = squeeze(sum(L_Xi(:, agent_idx, :).^2, 1)); 
                    
                    % 计算论文中的控制系数
                    N_i = neighbor_count_per_agent(agent_idx);
                    et_coeff = (et_sigma * et_a * (1 - et_a * N_i) / N_i);
                    
                    % 针对“每个诱导点”生成独立的阈值向量
                    Threshold = et_coeff * norm_z_sq + Bound;
                    
                    % 找出当前 agent 需要触发的诱导点索引
                    trigger_idx = find(norm_e_sq > Threshold);
                    
                    if dac_iter == 1
                        trigger_idx = 1:num_inducing_points; % 第一步全体广播初始化
                    end
                    
                    if ~isempty(trigger_idx)
                        % 仅更新“超标”的那些诱导点的广播状态
                        Xi_last_trigger(:, agent_idx, trigger_idx) = Xi(:, agent_idx, trigger_idx);
                        % 独立通信次数 +1
                        trigger_count_per_agent(agent_idx, trigger_idx) = trigger_count_per_agent(agent_idx, trigger_idx) + 1;
                    end
                end

                %% Step 4：系统整体收敛判断
                if max(abs(Xi(:) - Xi_prev(:))) < 1e-5
                    break;
                end
            end

            iter_converge = dac_iter;
            % 计算物理意义上的平均通信次数：平均“每个诱导点”被广播了多少次
            comm_train    = mean(trigger_count_per_agent(:));
            comm_test     = 0;
            fprintf('  收敛步数: %d  平均单点通信次数: %.1f 次\n', iter_converge, comm_train);

            % =======================================================
            % 核心修复 4: 提取最终聚合结果
            % 根据标准 DAC 定理，状态 Xi 会直接收敛到初始状态 Pi 的平均值！
            % 因此直接使用 Xi 作为共识结果，无需再做减法
            % =======================================================
            Xi_consensus = Xi;  
            phi_aggregated = zeros(output_dim, num_inducing_points);
            for dim_idx = 1:output_dim
                % 取第一个 agent 的收敛值即可（共识后各个智能体的值趋于一致）
                Xi_numerator   = squeeze(Xi_consensus(2*dim_idx-1, 1, :))';
                Xi_denominator = squeeze(Xi_consensus(2*dim_idx,   1, :))';
                
                if ismember(lower(current_method), {'gpoe','poe','bcm','rbcm'})
                    phi_aggregated(dim_idx,:) = Xi_numerator ./ max(Xi_denominator, eps);
                else
                    phi_aggregated(dim_idx,:) = Xi_numerator / num_agents;
                end
            end

        case ac_methods
            %% --- IP-AC：诱导点平均共识（一轮通信，无迭代）---
            comm_train    = 1;
            comm_test     = 0;
            iter_converge = 1;
            base_method   = strrep(lower(current_method), '_ac', '');

            phi_aggregated = zeros(output_dim, num_inducing_points);
            for inducing_idx = 1:num_inducing_points
                for dim_idx = 1:output_dim
                    mu_all_agents  = squeeze(local_mu_at_inducing(:, dim_idx, inducing_idx));
                    var_all_agents = squeeze(local_var_at_inducing(:, dim_idx, inducing_idx));
                    beta_all       = max(0.5*(log(prior_var)-log(var_all_agents)), eps);
                    switch base_method
                        case 'moe'
                            phi_aggregated(dim_idx, inducing_idx) = mean(mu_all_agents);
                        case 'gpoe'
                            prec = sum(beta_all ./ var_all_agents);
                            phi_aggregated(dim_idx, inducing_idx) = sum(beta_all .* mu_all_agents ./ var_all_agents) / max(prec, eps);
                        case 'poe'
                            prec = sum(1 ./ var_all_agents);
                            phi_aggregated(dim_idx, inducing_idx) = sum(mu_all_agents ./ var_all_agents) / max(prec, eps);
                        case 'bcm'
                            prec = sum(1 ./ var_all_agents) - (num_agents-1) / prior_var;
                            phi_aggregated(dim_idx, inducing_idx) = sum(mu_all_agents ./ var_all_agents) / max(prec, eps);
                        case 'rbcm'
                            prec = sum(beta_all ./ var_all_agents) + (1 - sum(beta_all)) / prior_var;
                            phi_aggregated(dim_idx, inducing_idx) = sum(beta_all .* mu_all_agents ./ var_all_agents) / max(prec, eps);
                    end
                end
            end
    end

    %% 8. 用聚合后的诱导点值训练 MaskedGP，然后批量预测
    MaskedGP = LocalGP_MultiOutput(input_dim, output_dim, num_inducing_points, 1e-6, SigmaF, SigmaL);
    MaskedGP.add_Alldata(InducingPoints_X, phi_aggregated);
    t_train_total = t_train_gp + toc;

    % 批量推断（利用已训练好的 Cholesky 分解，避免重复计算）
    tic;
    num_inducing_used = MaskedGP.DataQuantity;
    alpha_vec         = MaskedGP.alpha(1:num_inducing_used, :);       % [M × output_dim]
    chol_lower        = MaskedGP.L(1:num_inducing_used, 1:num_inducing_used); % [M × M] 下三角
    K_test_inducing   = MaskedGP.kernel(MaskedGP.X(:,1:num_inducing_used), X_eval'); % [M × num_eval]
    mu_normalized     = (alpha_vec' * K_test_inducing)';              % [num_eval × output_dim]
    V_chol_solve      = chol_lower \ K_test_inducing;                 % [M × num_eval]
    var_normalized    = max(SigmaF^2 - sum(V_chol_solve.^2, 1)', SigmaN^2); % [num_eval × 1]
    mu_pred           = mu_normalized .* repmat(Y_std, num_eval, 1) + repmat(Y_mean, num_eval, 1);
    var_pred          = repmat(var_normalized, 1, output_dim) .* repmat(Y_std.^2, num_eval, 1);
    t_test_total = toc;

    %% 9. 误差计算与保存
    prediction_error  = Y_eval - mu_pred;
    smse = mean(mean(prediction_error.^2) ./ Y_var_base);
    rmse = mean(sqrt(mean(prediction_error.^2)));
    nlpd = mean(mean(0.5*(log(2*pi*var_pred) + prediction_error.^2 ./ var_pred)));

    t_train_per_point = (t_train_total / num_train) * 1000;
    t_test_per_point  = (t_test_total  / num_eval)  * 1000;

    fprintf('  SMSE=%.4f  RMSE=%.4f  NLPD=%.4f  Train:%.2fms/pt  Test:%.2fms/pt\n', ...
        smse, rmse, nlpd, t_train_per_point, t_test_per_point);

    err_sq_mean = mean(prediction_error.^2, 2);
    smse_curve  = cumsum(err_sq_mean) ./ (1:num_eval)' / mean(Y_var_base);
    rmse_curve  = sqrt(cumsum(err_sq_mean) ./ (1:num_eval)');

    save(fullfile(SaveFolder, sprintf('%s_tr%d_mc%d.mat', current_method, tr_tag, seed)), ...
        'smse', 'rmse', 'nlpd', 't_train_total', 't_test_total', ...
        't_train_per_point', 't_test_per_point', ...
        'comm_train', 'comm_test', 'iter_converge', ...
        'current_method', 'seed', 'train_ratio', 'smse_curve', 'rmse_curve');
end
end