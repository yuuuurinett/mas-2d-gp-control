function run_testpoint_dataset(DatasetName, CurrentMode, train_ratio, seed, NumInducingPoints)
%RUN_TESTPOINT_DATASET Fast vectorized TP-DAC / TP-AC dataset evaluation.
%
% Important logic:
%   TP-DAC:
%       Pi is the test-point local reference signal.
%       Zeta is the DAC auxiliary state.
%       Xi = Pi - Zeta is the local consensus state.
%
%   TP-AC:
%       Pi is directly used as the initial static consensus state:
%           Xi_i(0) = Pi_i
%       The AC iteration converges to the initial average:
%           Xi_i(k) -> mean_i Xi_i(0) = mean_i Pi_i
%       No Zeta is used in TP-AC.
%
% The implementation vectorizes:
%   1) construction of P_info_matrix,
%   2) Laplacian multiplication over the agent dimension,
%   3) extraction of aggregated mean and variance.

if nargin < 3, train_ratio = 0.4; end
if nargin < 4, seed = 1; end
if nargin < 5, NumInducingPoints = 2500; end  % default matches max M
rng(seed);

fprintf('\n[测试点] %s  seed=%d  tr=%.0f%%\n', DatasetName, seed, train_ratio*100);

%% 1. 加载数据集
switch upper(DatasetName)
    case 'KIN40K'
        tr = load('KIN40K_train.mat');
        te = load('KIN40K_test.mat');
        hp = load('KIN40K_Hyperparameter.mat');
        train_x = tr.x; train_y = tr.y;
        test_x  = te.xtest; test_y = te.ytest;

    case 'POL'
        tr = load(fullfile('POL','POL_train.mat'));
        te = load(fullfile('POL','POL_test.mat'));
        hp = load(fullfile('POL','POL_Hyperparameter.mat'));
        train_x = tr.x; train_y = tr.y;
        test_x  = te.xtest; test_y = te.ytest;

    case 'PUMADYN32NM'
        tr = load(fullfile('PUMADYN32NM','PUMADYN32NM_train.mat'));
        te = load(fullfile('PUMADYN32NM','PUMADYN32NM_test.mat'));
        hp = load(fullfile('PUMADYN32NM','PUMADYN32NM_Hyperparameter.mat'));
        train_x = tr.x; train_y = tr.y;
        test_x  = te.xtest; test_y = te.ytest;

    case 'SARCOS'
        tr = load(fullfile('SARCOS','SARCOS_train.mat'));
        te = load(fullfile('SARCOS','SARCOS_test.mat'));
        hp_raw = load(fullfile('SARCOS','SARCOS_GP_Hyperparameter.mat'));

        hp.SigmaF = mean(cell2mat(hp_raw.SigmaF_set));
        hp.SigmaN = mean(cell2mat(hp_raw.SigmaN_set));
        hp.SigmaL = mean(cell2mat(hp_raw.SigmaL_set'), 2);

        train_x = tr.sarcos_inv(:,1:21);
        train_y = tr.sarcos_inv(:,22:28);
        test_x  = te.sarcos_inv_test(:,1:21);
        test_y  = te.sarcos_inv_test(:,22:28);

    otherwise
        error('未知数据集: %s', DatasetName);
end

%% 2. 归一化
if size(hp.SigmaL,1) > 1 && size(hp.SigmaL,2) > 1
    hp.SigmaL = mean(hp.SigmaL, 1);
end
SigmaL = hp.SigmaL(:);

if numel(hp.SigmaF) > 1, hp.SigmaF = mean(hp.SigmaF); end
if numel(hp.SigmaN) > 1, hp.SigmaN = mean(hp.SigmaN); end
SigmaF = hp.SigmaF;
SigmaN = hp.SigmaN;

num_train_samples = round(size(train_x,1) * train_ratio);
train_indices = randperm(size(train_x,1), num_train_samples);
X_train = train_x(train_indices,:);
Y_train = train_y(train_indices,:);

test_indices = randperm(size(test_x,1));
X_test = test_x(test_indices,:);
Y_test = test_y(test_indices,:);

X_mean = mean(X_train,1);
X_std  = std(X_train,0,1);
X_std(X_std == 0) = 1;

if ~(max(abs(X_mean)) < 1e-2 && max(abs(X_std - 1)) < 1e-2)
    X_train = (X_train - X_mean) ./ X_std;
    X_test  = (X_test  - X_mean) ./ X_std;
    SigmaL  = SigmaL ./ X_std(:);
end

Y_mean = mean(Y_train,1);
Y_std  = std(Y_train,0,1);
Y_std(Y_std == 0) = 1;

if max(abs(Y_mean)) < 1e-2 && max(abs(Y_std - 1)) < 1e-2
    Y_mean = zeros(1,size(Y_train,2));
    Y_std  = ones(1,size(Y_train,2));
else
    Y_train = (Y_train - Y_mean) ./ Y_std;
    SigmaF = SigmaF / mean(Y_std);
    SigmaN = SigmaN / mean(Y_std);
end

prior_var = SigmaF^2;

[num_train, input_dim] = size(X_train);
output_dim = size(Y_train,2);
num_eval = min(NumInducingPoints, size(X_test,1));  % match inducing point count M

X_eval = X_test(1:num_eval,:);
Y_eval = Y_test(1:num_eval,:);
Y_var_base = var(Y_eval,0,1);
Y_var_base_safe = Y_var_base;
Y_var_base_safe(Y_var_base_safe <= eps) = eps;

fprintf('Train=%d  Test=%d  x_dim=%d  y_dim=%d\n', ...
    num_train, num_eval, input_dim, output_dim);

%% 3. 分布式系统参数
num_agents = 6;
dac_step_size = 0.01;
dac_gain = 10;
max_iters = 3000;
conv_tol = 1e-5;

max_data_per_agent = min(floor(num_train / num_agents), 3000);

MultiAgentSystem = Manipulator_2D_2DoF_SetMASTopology(num_agents, 1);
Laplacian = MultiAgentSystem.Agent_Topology.LaplacianMatrix;

neighbor_count_per_agent = sum(Laplacian < 0, 2);
num_directed_neighbor_links = sum(neighbor_count_per_agent);
neighbor_mask = (Laplacian < 0) | (Laplacian' < 0);
neighbor_mask(1:size(neighbor_mask,1)+1:end) = false;

%% 4. 训练局部 GP
tic_local_gp = tic;
local_gp_set = cell(num_agents,1);

for agent_idx = 1:num_agents
    data_idx = (agent_idx-1)*max_data_per_agent + 1 : ...
        min(agent_idx*max_data_per_agent, num_train);

    local_gp_set{agent_idx} = LocalGP_MultiOutput(input_dim, output_dim, ...
        max_data_per_agent, SigmaN, SigmaF, SigmaL);

    local_gp_set{agent_idx}.add_Alldata(X_train(data_idx,:)', Y_train(data_idx,:)');
    local_gp_set{agent_idx}.tau = 1e-8;
    local_gp_set{agent_idx}.delta = 0.01;
end

t_train_gp = toc(tic_local_gp);
t_tp_local_gp_train = t_train_gp;
fprintf('局部GP训练完成: %.4fs\n', t_train_gp);

%% 5. 方法列表
dac_methods = {'poe','gpoe','moe','bcm','rbcm'};
ac_methods  = {'moe_ac','gpoe_ac','poe_ac','bcm_ac','rbcm_ac'};

if strcmpi(CurrentMode,'all')
    AllModes = [dac_methods, ac_methods];
else
    AllModes = {lower(CurrentMode)};
end

SaveFolder = fullfile('Result','Dataset',DatasetName);
if ~exist(SaveFolder,'dir'), mkdir(SaveFolder); end

p_dim = 2 * output_dim;
tr_tag = round(train_ratio * 100);

%% 6. 批量预计算测试点上的局部 GP 预测
fprintf('[预计算] 正在批量计算 %d 个测试点的局部预测...\n', num_eval);
tic_test_local_prediction = tic;

local_mu_at_testpoints  = zeros(num_agents, output_dim, num_eval);
local_var_at_testpoints = zeros(num_agents, output_dim, num_eval);

for agent_idx = 1:num_agents
    % batch_predict_external should return:
    %   mu_batch  [num_eval x output_dim]
    %   var_batch [num_eval x 1] or [num_eval x output_dim]
    [mu_batch, var_batch] = batch_predict_external( ...
        local_gp_set{agent_idx}, X_eval', SigmaN, SigmaF);

    local_mu_at_testpoints(agent_idx,:,:) = permute(mu_batch, [3 2 1]);

    if isvector(var_batch)
        local_var_at_testpoints(agent_idx,:,:) = ...
            permute(repmat(var_batch(:), 1, output_dim), [3 2 1]);
    else
        local_var_at_testpoints(agent_idx,:,:) = permute(var_batch, [3 2 1]);
    end
end

precompute_time = toc(tic_test_local_prediction);
fprintf('预计算完成: %.2fs\n', precompute_time);

%% 7. 主循环
for method_idx = 1:numel(AllModes)
    current_method = AllModes{method_idx};
    base_method_name = strrep(lower(current_method), '_ac', '');

    fprintf('\n[%d/%d] %s\n', method_idx, numel(AllModes), current_method);
    t_method_start = tic;

    % When running all methods in one call, local test predictions are
    % precomputed once and reused. We allocate that shared cost evenly
    % across methods for per-method timing breakdown. If only one method
    % is run, this equals the full precomputation time.
    t_tp_test_local_prediction = precompute_time / max(numel(AllModes), 1);

    %% 7.1 构建 P 信息矩阵，维度 [p_dim x num_agents x num_eval]
    tic_build_info = tic;
    P_info_matrix = build_information_matrix( ...
        base_method_name, ...
        local_mu_at_testpoints, ...
        local_var_at_testpoints, ...
        SigmaN, prior_var, num_agents, output_dim);
    t_tp_build_information = toc(tic_build_info);

    %% 7.2 TP-DAC 或 TP-AC 共识
    if ismember(lower(current_method), dac_methods)
        tic_consensus = tic;
        % ------------------------------------------------------------
        % TP-DAC:
        %   Xi = Pi - Zeta.
        %   Zeta_dot = kappa * L * Xi.
        %   Here Pi is fixed over the test-point consensus iteration,
        %   but we keep the DAC structure to match the TP-DAC method.
        % ------------------------------------------------------------
        Zeta = zeros(size(P_info_matrix));
        iter_converge = max_iters;

        for iter = 1:max_iters
            Zeta_prev = Zeta;

            Xi_now = P_info_matrix - Zeta;
            L_Xi = laplacian_multiply_agent_dim(Xi_now, Laplacian);

            Zeta = Zeta + dac_step_size * dac_gain * L_Xi;

            if max(abs(Zeta(:) - Zeta_prev(:))) < conv_tol
                iter_converge = iter;
                break;
            end
        end

        Xi_consensus = P_info_matrix - Zeta;
        comm_test = iter_converge;
        comm_train = 0;

        t_tp_consensus = toc(tic_consensus);
        fprintf('  [TP-DAC] 收敛步数:%d  通信轮数:%d\n', iter_converge, comm_test);

    else
        tic_consensus = tic;
        % ------------------------------------------------------------
        % TP-AC:
        %   Xi_i(0) = Pi_i.
        %   Xi_dot = -kappa * L * Xi.
        %   Therefore Xi_i converges to the initial average mean_i Pi_i.
        %
        % Important:
        %   AC must NOT use closed-form aggregation here if we want to
        %   evaluate test-point AC communication/convergence. It must
        %   iterate from the initial local states Pi to their average.
        % ------------------------------------------------------------
        Xi = P_info_matrix;
        Xi_initial_average = mean(P_info_matrix, 2); %#ok<NASGU>
        iter_converge = max_iters;

        for iter = 1:max_iters
            Xi_prev = Xi;

            L_Xi = laplacian_multiply_agent_dim(Xi, Laplacian);
            Xi = Xi - dac_step_size * dac_gain * L_Xi;

            if max(abs(Xi(:) - Xi_prev(:))) < conv_tol
                iter_converge = iter;
                break;
            end
        end

        Xi_consensus = Xi;
        comm_test = iter_converge;
        comm_train = 0;

        t_tp_consensus = toc(tic_consensus);
        fprintf('  [TP-AC] 收敛步数:%d  通信轮数:%d\n', iter_converge, comm_test);
    end

    %% 7.3 从最终 consensus state 中提取聚合均值和方差
    tic_post_consensus = tic;
    final_mu_pred = extract_mean_from_consensus( ...
        base_method_name, Xi_consensus, output_dim, num_agents);

    %% 7.4 向量化计算聚合方差
    final_var_pred = aggregate_variance_vectorized( ...
        base_method_name, local_var_at_testpoints, ...
        SigmaN, prior_var, num_agents, output_dim);

    % Time for extracting the final prediction from the consensus state.
    t_tp_post_consensus = toc(tic_post_consensus);

    %% 7.5 反归一化、误差计算与保存
    % End-to-end test time for this TP method. This includes the allocated
    % local test prediction precomputation time and the method-specific work
    % after the loop starts.
    t_test_total = precompute_time / max(numel(AllModes),1) + toc(t_method_start);

    % Explicit timing breakdown used for later analysis. The small difference
    % between t_test_total and the sum below is recorded as method overhead.
    t_tp_total_test = t_tp_test_local_prediction + ...
                      t_tp_build_information + ...
                      t_tp_consensus + ...
                      t_tp_post_consensus;

    t_tp_method_overhead = t_test_total - t_tp_total_test;

    mu_pred = final_mu_pred .* repmat(Y_std, num_eval, 1) + repmat(Y_mean, num_eval, 1);
    var_pred = final_var_pred .* repmat(Y_std.^2, num_eval, 1);
    var_pred = max(var_pred, eps);

    prediction_error = Y_eval - mu_pred;

    smse = mean(mean(prediction_error.^2) ./ Y_var_base_safe);
    rmse = mean(sqrt(mean(prediction_error.^2)));
    nlpd = mean(mean(0.5 * (log(2*pi*var_pred) + prediction_error.^2 ./ var_pred)));

    t_train_per_point = (t_train_gp / num_train) * 1000;
    t_test_per_point  = (t_test_total / num_eval) * 1000;

    fprintf('  SMSE=%.4f  RMSE=%.4f  NLPD=%.4f  Train:%.2fms/pt  Test:%.2fms/pt\n', ...
        smse, rmse, nlpd, t_train_per_point, t_test_per_point);

    err_sq_mean = mean(prediction_error.^2, 2);
    smse_curve = cumsum(err_sq_mean) ./ (1:num_eval)' / mean(Y_var_base_safe);
    rmse_curve = sqrt(cumsum(err_sq_mean) ./ (1:num_eval)');

    event_count_test = comm_test;
    trigger_ratio_test = comm_test / max(iter_converge, eps);

    fprintf('  Timing TP: LocalPred=%.3fs  BuildInfo=%.3fs  Consensus=%.3fs  Post=%.3fs  Total=%.3fs\n', ...
        t_tp_test_local_prediction, t_tp_build_information, t_tp_consensus, t_tp_post_consensus, t_tp_total_test);

    trigger_count_per_agent = zeros(num_agents, num_eval);

    save(fullfile(SaveFolder, sprintf('%s_tp_tr%d_mc%d.mat', current_method, tr_tag, seed)), ...
        'smse', 'rmse', 'nlpd', ...
        't_train_gp', 't_test_total', ...
        't_train_per_point', 't_test_per_point', ...
        't_tp_local_gp_train', 't_tp_test_local_prediction', ...
        't_tp_build_information', 't_tp_consensus', ...
        't_tp_post_consensus', 't_tp_method_overhead', ...
        't_tp_total_test', ...
        'comm_train', 'comm_test', 'iter_converge', ...
        'neighbor_count_per_agent', 'num_directed_neighbor_links', 'neighbor_mask', ...
        'trigger_count_per_agent', 'event_count_test', 'trigger_ratio_test', ...
        'current_method', 'seed', 'train_ratio', ...
        'smse_curve', 'rmse_curve');
end

fprintf('\n[%s] test-point evaluation done. seed=%d tr=%d%%\n\n', ...
    DatasetName, seed, tr_tag);

end

%% ========================================================================
% Helper functions
% ========================================================================

function P_info_matrix = build_information_matrix(method_name, Mu, Var, SigmaN, prior_var, num_agents, output_dim)
% Mu, Var: [num_agents x output_dim x num_eval]
% P_info_matrix: [2*output_dim x num_agents x num_eval]

[~, ~, num_eval] = size(Mu);
p_dim = 2 * output_dim;
P_info_matrix = zeros(p_dim, num_agents, num_eval);

for dim_idx = 1:output_dim
    mu_d  = squeeze(Mu(:, dim_idx, :));   % [num_agents x num_eval]
    var_d = squeeze(Var(:, dim_idx, :));  % [num_agents x num_eval]

    var_d = max(var_d, SigmaN^2);
    beta_d = max(min(0.5 * (log(prior_var) - log(var_d)), 10), eps);

    switch lower(method_name)
        case 'moe'
            numerator   = num_agents * mu_d;
            denominator = num_agents * (var_d + mu_d.^2);

        case 'gpoe'
            numerator   = num_agents * beta_d .* mu_d ./ var_d;
            denominator = num_agents * beta_d ./ var_d;

        case 'poe'
            numerator   = num_agents * mu_d ./ var_d;
            denominator = num_agents ./ var_d;

        case 'bcm'
            numerator   = num_agents * mu_d ./ var_d;
            denominator = num_agents ./ var_d - (num_agents - 1) / prior_var;

        case 'rbcm'
            numerator   = num_agents * beta_d .* mu_d ./ var_d;
            denominator = num_agents * beta_d ./ var_d + (1 - num_agents * beta_d) / prior_var;

        otherwise
            error('Unknown aggregation method: %s', method_name);
    end

    P_info_matrix(2*dim_idx-1,:,:) = reshape(numerator,   [1, num_agents, num_eval]);
    P_info_matrix(2*dim_idx,  :,:) = reshape(denominator, [1, num_agents, num_eval]);
end
end

function L_X = laplacian_multiply_agent_dim(X, L)
% Applies the graph Laplacian over the agent dimension only.
%
% X has size [p_dim x num_agents x num_points].
% For each information dimension r and point m, this computes:
%   L_X(r,:,m) = X(r,:,m) * L'
% Therefore:
%   L_X(:,i,m) = sum_j L(i,j) * X(:,j,m)
% which equals (L * x)_i over the agent dimension.

[p_dim, num_agents, num_points] = size(X);

X_flat = reshape(permute(X, [1 3 2]), p_dim*num_points, num_agents);
LX_flat = X_flat * L';

L_X = ipermute(reshape(LX_flat, p_dim, num_points, num_agents), [1 3 2]);
end

function final_mu_pred = extract_mean_from_consensus(method_name, Xi_consensus, output_dim, num_agents)
% Xi_consensus: [p_dim x num_agents x num_eval]
% Output: [num_eval x output_dim]

num_eval = size(Xi_consensus, 3);
final_mu_pred = zeros(num_eval, output_dim);

Xi_mean = squeeze(mean(Xi_consensus, 2)); % [p_dim x num_eval]
if num_eval == 1
    Xi_mean = reshape(Xi_mean, [], 1);
end

for dim_idx = 1:output_dim
    xi1 = Xi_mean(2*dim_idx-1, :)';
    xi2 = Xi_mean(2*dim_idx,   :)';

    if ismember(lower(method_name), {'gpoe','poe','bcm','rbcm'})
        final_mu_pred(:, dim_idx) = xi1 ./ max(xi2, eps);
    else
        % MoE:
        % P_info numerator was num_agents * mu_i.
        % After average consensus, xi1 = sum_i mu_i.
        % Therefore divide by num_agents to obtain mean_i mu_i.
        final_mu_pred(:, dim_idx) = xi1 / num_agents;
    end
end
end

function final_var_pred = aggregate_variance_vectorized(method_name, Var, SigmaN, prior_var, num_agents, output_dim)
% Var: [num_agents x output_dim x num_eval]
% Output: [num_eval x output_dim]

num_eval = size(Var, 3);
final_var_pred = zeros(num_eval, output_dim);

for dim_idx = 1:output_dim
    var_d = squeeze(Var(:, dim_idx, :)); % [num_agents x num_eval]
    var_d = max(var_d, SigmaN^2);

    beta_d = max(0.5 * (log(prior_var) - log(var_d)), eps);

    switch lower(method_name)
        case 'moe'
            var_out = mean(var_d, 1);

        case 'gpoe'
            precision = sum(beta_d ./ var_d, 1);
            var_out = 1 ./ max(precision, eps);

        case 'poe'
            precision = sum(1 ./ var_d, 1);
            var_out = 1 ./ max(precision, eps);

        case 'bcm'
            precision = sum(1 ./ var_d, 1) - (num_agents - 1) / prior_var;
            var_out = 1 ./ max(precision, eps);

        case 'rbcm'
            precision = sum(beta_d ./ var_d, 1) + (1 - sum(beta_d, 1)) / prior_var;
            var_out = 1 ./ max(precision, eps);

        otherwise
            error('Unknown aggregation method: %s', method_name);
    end

    final_var_pred(:, dim_idx) = max(var_out(:), SigmaN^2);
end
end