function run_testpoint_dataset(DatasetName, CurrentMode, train_ratio, seed)

if nargin < 3, train_ratio = 0.4; end
if nargin < 4, seed = 1;          end
rng(seed);

fprintf('\n[测试点] %s  seed=%d  tr=%.0f%%\n', DatasetName, seed, train_ratio*100);

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
    otherwise, error('未知数据集: %s', DatasetName);
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
output_dim = size(Y_train, 2);
num_eval   = min(3000, size(X_test,1));
X_eval     = X_test(1:num_eval,:);
Y_eval     = Y_test(1:num_eval,:);
Y_var_base = var(Y_eval, 0, 1);
fprintf('Train=%d  Test=%d  x_dim=%d  y_dim=%d\n', num_train, num_eval, input_dim, output_dim);

%% 3. 分布式系统参数
num_agents         = 6;
dac_step_size      = 0.01;
dac_gain           = 10;
max_data_per_agent = min(floor(num_train / num_agents), 3000);

MultiAgentSystem = Manipulator_2D_2DoF_SetMASTopology(num_agents, 1);
Laplacian        = MultiAgentSystem.Agent_Topology.LaplacianMatrix;

neighbor_count_per_agent = sum(Laplacian < 0, 2);
max_neighbor_count       = max(neighbor_count_per_agent);

% ET 参数
et_sigma = 0.2;
et_a     = 0.9 / max_neighbor_count;
fprintf('ET 参数: sigma=%.2f  a=%.4f  (max|N_i|=%d)\n', et_sigma, et_a, max_neighbor_count);

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

%% 5. 方法列表
dac_methods = {'moe','gpoe','poe','bcm','rbcm'};
ac_methods  = {'moe_ac','gpoe_ac','poe_ac','bcm_ac','rbcm_ac'};
if strcmpi(CurrentMode,'all')
    AllModes = [dac_methods, ac_methods];
else
    AllModes = {lower(CurrentMode)};
end

SaveFolder = fullfile('Result','Dataset',DatasetName);
if ~exist(SaveFolder,'dir'), mkdir(SaveFolder); end
p_dim  = 2 * output_dim;
tr_tag = round(train_ratio * 100);

%% 6. 批量预计算：对所有测试点做局部 GP 预测
% local_mu_at_testpoints  [num_agents × output_dim × num_eval]：各 agent 在各测试点的预测均值
% local_var_at_testpoints [num_agents × output_dim × num_eval]：各 agent 在各测试点的预测方差
fprintf('\n[批量预计算] 正在计算 %d 个测试点的局部预测...\n', num_eval);
tic;
local_mu_at_testpoints  = zeros(num_agents, output_dim, num_eval);
local_var_at_testpoints = zeros(num_agents, output_dim, num_eval);
X_eval_col = X_eval';  % 转为列格式 [input_dim × num_eval]

for agent_idx = 1:num_agents
    [mu_batch, var_batch] = batch_predict_external( ...
        local_gp_set{agent_idx}, X_eval_col, SigmaN, SigmaF);
    for test_idx = 1:num_eval
        local_mu_at_testpoints(agent_idx,:,test_idx)  = mu_batch(test_idx,:);
        local_var_at_testpoints(agent_idx,:,test_idx) = repmat(var_batch(test_idx), 1, output_dim);
    end
end
precompute_time = toc;
fprintf('预计算完成，耗时: %.2fs (%.4f ms/pt)\n', precompute_time, precompute_time/num_eval*1000);

%% 7. 聚合方法主循环
for method_idx = 1:numel(AllModes)
    current_method   = AllModes{method_idx};
    base_method_name = strrep(lower(current_method), '_ac', '');
    fprintf('\n[%d/%d] 方法: %s\n', method_idx, numel(AllModes), current_method);

    final_mu_pred  = zeros(num_eval, output_dim);
    final_var_pred = zeros(num_eval, output_dim);
    tic;

    if ismember(lower(current_method), dac_methods)
        %% --- TP-DAC：测试点分布式平均共识（DAC）---

        % 构建信息向量矩阵 P_info_matrix
        % P_info_matrix [p_dim × num_agents × num_eval]：每个测试点上每个 agent 的信息向量
        P_info_matrix = zeros(p_dim, num_agents, num_eval);
        for test_idx = 1:num_eval
            for agent_idx = 1:num_agents
                local_mu_i  = local_mu_at_testpoints(agent_idx,:,test_idx)';
                local_var_i = local_var_at_testpoints(agent_idx,:,test_idx)';
                for dim_idx = 1:output_dim
                    safe_var  = max(local_var_i(dim_idx), SigmaN^2);
                    beta_gpoe = max(min(0.5*(log(prior_var)-log(safe_var)), 10), eps);
                    switch lower(current_method)
                        case 'moe'
                            P_info_matrix(2*dim_idx-1, agent_idx, test_idx) = num_agents * local_mu_i(dim_idx);
                            P_info_matrix(2*dim_idx,   agent_idx, test_idx) = num_agents * (safe_var + local_mu_i(dim_idx)^2);
                        case 'gpoe'
                            P_info_matrix(2*dim_idx-1, agent_idx, test_idx) = num_agents * beta_gpoe * local_mu_i(dim_idx) / safe_var;
                            P_info_matrix(2*dim_idx,   agent_idx, test_idx) = num_agents * beta_gpoe / safe_var;
                        case 'poe'
                            P_info_matrix(2*dim_idx-1, agent_idx, test_idx) = num_agents * local_mu_i(dim_idx) / safe_var;
                            P_info_matrix(2*dim_idx,   agent_idx, test_idx) = num_agents / safe_var;
                        case 'bcm'
                            P_info_matrix(2*dim_idx-1, agent_idx, test_idx) = num_agents * local_mu_i(dim_idx) / safe_var;
                            P_info_matrix(2*dim_idx,   agent_idx, test_idx) = num_agents / safe_var - (num_agents-1) / prior_var;
                        case 'rbcm'
                            P_info_matrix(2*dim_idx-1, agent_idx, test_idx) = num_agents * beta_gpoe * local_mu_i(dim_idx) / safe_var;
                            P_info_matrix(2*dim_idx,   agent_idx, test_idx) = num_agents * beta_gpoe / safe_var + (1 - num_agents*beta_gpoe) / prior_var;
                    end
                end
            end
        end

        % DAC 迭代
       % =========================================================================
        % [核心修复] TP-DAC 纯正分布式平均共识与独立事件触发
        % =========================================================================
        Xi              = P_info_matrix;  % 真实连续状态 [p_dim, num_agents, num_eval]
        Xi_last_trigger = P_info_matrix;  % 广播状态快照
        
        % 独立通信计数器：记录 [每个智能体] 对 [每个测试点] 的通信触发次数
        trigger_count_per_agent = zeros(num_agents, num_eval); 
        
        dac_iter  = 0;
        max_iters = 3000;
        
        % 防 Zeno 现象的衰减 Bound 参数
        c_0 = 0.1;           
        alpha_decay = 0.02;  
        comm_train = 0;

        while dac_iter < max_iters
            dac_iter  = dac_iter + 1;
            Xi_prev = Xi;

            %% Step 1：计算邻居的广播状态驱动力 (L * Xi_last_trigger)
            L_Xi = zeros(size(Xi_last_trigger));
            for agent_idx = 1:num_agents
                L_Xi(:, agent_idx, :) = sum(Xi_last_trigger .* reshape(Laplacian(agent_idx,:), 1, num_agents, 1), 2);
            end

            %% Step 2：真实状态演化 (\dot{Xi} = - L * \hat{Xi})
            for agent_idx = 1:num_agents
                Xi(:, agent_idx, :) = Xi(:, agent_idx, :) - dac_step_size * dac_gain * L_Xi(:, agent_idx, :);
            end

            %% Step 3：事件触发判断 (Point-wise 完全解耦，各测试点独立触发)
            Bound = c_0 * exp(-alpha_decay * dac_iter); 
            
            for agent_idx = 1:num_agents
                % 误差矩阵 [p_dim, 1, num_eval]
                error_i = Xi_last_trigger(:, agent_idx, :) - Xi(:, agent_idx, :);
                
                % 沿着特征维度求平方和，保留每个测试点的独立性 (降维为 num_eval 长度的向量)
                norm_e_sq = squeeze(sum(error_i.^2, 1));    
                norm_z_sq = squeeze(sum(L_Xi(:, agent_idx, :).^2, 1)); 
                
                % 统一转为列向量防止维度不匹配
                if size(norm_e_sq, 1) == 1, norm_e_sq = norm_e_sq'; end
                if size(norm_z_sq, 1) == 1, norm_z_sq = norm_z_sq'; end
                
                % 计算当前 agent 的控制系数
                N_i = neighbor_count_per_agent(agent_idx);
                et_coeff = (et_sigma * et_a * (1 - et_a * N_i) / N_i);
                
                % 独立测试点的阈值向量
                Threshold = et_coeff * norm_z_sq + Bound;
                
                % 找出当前 agent 需要触发的“测试点”索引
                trigger_idx = find(norm_e_sq > Threshold);
                
                if dac_iter == 1
                    trigger_idx = 1:num_eval; % 第一步强制全体触发广播以初始化
                end
                
                if ~isempty(trigger_idx)
                    % 仅更新触发了的测试点状态
                    Xi_last_trigger(:, agent_idx, trigger_idx) = Xi(:, agent_idx, trigger_idx);
                    trigger_count_per_agent(agent_idx, trigger_idx) = trigger_count_per_agent(agent_idx, trigger_idx) + 1;
                end
            end

            %% Step 4：收敛判断
            if max(abs(Xi(:) - Xi_prev(:))) < 1e-5
                break;
            end
        end

        iter_converge = dac_iter;
        comm_test     = mean(trigger_count_per_agent(:));
        fprintf('    -> 收敛步数: %d  平均ET单点触发次数: %.1f\n', iter_converge, comm_test);

        %% Step 5：提取聚合预测
        Xi_consensus = Xi;
        for test_idx = 1:num_eval
            for dim_idx = 1:output_dim
                Xi_num = Xi_consensus(2*dim_idx-1, 1, test_idx);
                Xi_den = Xi_consensus(2*dim_idx,   1, test_idx);
                if ismember(lower(current_method), {'gpoe','poe','bcm','rbcm'})
                    final_mu_pred(test_idx, dim_idx) = Xi_num / max(Xi_den, eps);
                else
                    final_mu_pred(test_idx, dim_idx) = Xi_num / num_agents;
                end
            end
        end
    else
        %% --- TP-AC：测试点平均共识（一轮通信，无迭代）---
        comm_train    = 0;
        comm_test     = 1;
        iter_converge = 1;
        for test_idx = 1:num_eval
            for dim_idx = 1:output_dim
                mu_all_agents  = local_mu_at_testpoints(:, dim_idx, test_idx);
                var_all_agents = local_var_at_testpoints(:, dim_idx, test_idx);
                beta_all       = max(0.5*(log(prior_var)-log(var_all_agents)), eps);
                switch base_method_name
                    case 'moe'
                        final_mu_pred(test_idx, dim_idx) = mean(mu_all_agents);
                    case 'gpoe'
                        prec = sum(beta_all ./ var_all_agents);
                        final_mu_pred(test_idx, dim_idx) = sum(beta_all .* mu_all_agents ./ var_all_agents) / max(prec, eps);
                    case 'poe'
                        prec = sum(1 ./ var_all_agents);
                        final_mu_pred(test_idx, dim_idx) = sum(mu_all_agents ./ var_all_agents) / max(prec, eps);
                    case 'bcm'
                        prec = sum(1 ./ var_all_agents) - (num_agents-1) / prior_var;
                        final_mu_pred(test_idx, dim_idx) = sum(mu_all_agents ./ var_all_agents) / max(prec, eps);
                    case 'rbcm'
                        prec = sum(beta_all ./ var_all_agents) + (1 - sum(beta_all)) / prior_var;
                        final_mu_pred(test_idx, dim_idx) = sum(beta_all .* mu_all_agents ./ var_all_agents) / max(prec, eps);
                end
            end
        end
    end

    %% 8. 方差聚合（DAC 和 AC 共用同样的方差公式）
    for test_idx = 1:num_eval
        for dim_idx = 1:output_dim
            var_all_agents = local_var_at_testpoints(:, dim_idx, test_idx);
            beta_all       = max(0.5*(log(prior_var)-log(var_all_agents)), eps);
            switch base_method_name
                case 'moe'
                    final_var_pred(test_idx, dim_idx) = mean(var_all_agents);
                case 'gpoe'
                    prec = sum(beta_all ./ var_all_agents);
                    final_var_pred(test_idx, dim_idx) = 1 / max(prec, eps);
                case 'poe'
                    prec = sum(1 ./ var_all_agents);
                    final_var_pred(test_idx, dim_idx) = 1 / max(prec, eps);
                case 'bcm'
                    prec = sum(1 ./ var_all_agents) - (num_agents-1) / prior_var;
                    final_var_pred(test_idx, dim_idx) = 1 / max(prec, eps);
                case 'rbcm'
                    prec = sum(beta_all ./ var_all_agents) + (1 - sum(beta_all)) / prior_var;
                    final_var_pred(test_idx, dim_idx) = 1 / max(prec, eps);
            end
            final_var_pred(test_idx, dim_idx) = max(final_var_pred(test_idx, dim_idx), SigmaN^2);
        end
    end

    t_test_total = precompute_time / numel(AllModes) + toc;

    %% 9. 反归一化与保存
    mu_pred  = final_mu_pred  .* repmat(Y_std, num_eval, 1) + repmat(Y_mean, num_eval, 1);
    var_pred = final_var_pred .* repmat(Y_std.^2, num_eval, 1);

    prediction_error = Y_eval - mu_pred;
    smse = mean(mean(prediction_error.^2) ./ Y_var_base);
    rmse = mean(sqrt(mean(prediction_error.^2)));
    nlpd = mean(mean(0.5*(log(2*pi*var_pred) + prediction_error.^2 ./ var_pred)));

    t_train_per_point = (t_train_gp  / num_train) * 1000;
    t_test_per_point  = (t_test_total / num_eval)  * 1000;

    fprintf('  SMSE=%.4f  RMSE=%.4f  NLPD=%.4f  Train:%.2fms/pt  Test:%.2fms/pt\n', ...
        smse, rmse, nlpd, t_train_per_point, t_test_per_point);

    err_sq_mean = mean(prediction_error.^2, 2);
    smse_curve  = cumsum(err_sq_mean) ./ (1:num_eval)' / mean(Y_var_base);
    rmse_curve  = sqrt(cumsum(err_sq_mean) ./ (1:num_eval)');

    save(fullfile(SaveFolder, sprintf('%s_tp_tr%d_mc%d.mat', current_method, tr_tag, seed)), ...
        'smse', 'rmse', 'nlpd', 't_train_gp', 't_test_total', ...
        't_train_per_point', 't_test_per_point', ...
        'comm_train', 'comm_test', 'iter_converge', ...
        'current_method', 'seed', 'train_ratio', 'smse_curve', 'rmse_curve');
end
end