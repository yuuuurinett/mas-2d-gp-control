function [MaskedGP, Zeta_vector, Zeta_last_trigger, trigger_count] = gp_masked_aggregation_update( ...
    P, Zeta_vector, L, ...
    Kappa_P, AgentQuantity, NumInducingPoints, TimeStep, ...
    InducingPoints_Coordinates, SigmaF, SigmaL, x_dim, method, p_dim, ...
    Zeta_last_trigger, et_sigma, et_a, neighbor_count)
% 新增参数：
%   Zeta_last_trigger [p_dim × AgentQuantity × M]：上次触发时的广播值
%   et_sigma, et_a：ET 参数
%   neighbor_count [AgentQuantity × 1]：每个 agent 的邻居数量
% 新增输出：
%   Zeta_last_trigger：更新后的广播值
%   trigger_count [AgentQuantity × 1]：本步各 agent 的触发次数（0或1）

method    = lower(method);
M         = NumInducingPoints;
y_dim     = 2;
prior_var = SigmaF^2;

trigger_count = zeros(AgentQuantity, 1);

%% Step 1: DAC 更新（Euler 离散化）+ ET 触发判定
if TimeStep > 0

    % --- Euler 更新：对所有诱导点同时做一步 ---
    % Zeta_vector: [p_dim × AgentQuantity × M]
    % dZeta/dt = Kappa * (P - Zeta) * L'
    % 用 last_trigger 值参与共识（ET 核心：用广播值而非实时值）
    for InducingPointIdx = 1:M
        P_m             = P(:, :, InducingPointIdx);             % [p_dim × AgentQuantity]
        Zeta_m          = Zeta_vector(:, :, InducingPointIdx);   % [p_dim × AgentQuantity]
        Zeta_trigger_m  = Zeta_last_trigger(:, :, InducingPointIdx); % [p_dim × AgentQuantity]

        % 共识输入用广播值
        consensus_input = P_m - Zeta_trigger_m;  % [p_dim × AgentQuantity]
        dZeta = Kappa_P * consensus_input * L';   % [p_dim × AgentQuantity]
        Zeta_vector(:, :, InducingPointIdx) = Zeta_m + TimeStep * dZeta;
    end

    % --- ET 触发判定：对每个 agent，合并所有诱导点 ---
    % 对应 Dimarogonas 2012 公式(11)：
    %   ||e_i||² ≤ σ·a·(1-a·|N_i|)/|N_i| · ||z_i||²
    for agent_idx = 1:AgentQuantity
        % e_i = Zeta_last_trigger_i - Zeta_i（上次广播值与当前值之差）
        % 展平所有诱导点：[p_dim × M] → [p_dim*M × 1]
        error_i = Zeta_last_trigger(:, agent_idx, :) - Zeta_vector(:, agent_idx, :);
        error_i_flat = reshape(error_i, [], 1);

        % z_i = sum_{j∈N_i}(Zeta_i - Zeta_j)，用实时值计算
        neighbor_mask = reshape(L(agent_idx,:) < 0, 1, AgentQuantity, 1);
        relative_diff = sum((Zeta_vector(:, agent_idx, :) - Zeta_vector) .* neighbor_mask, 2);
        relative_diff_flat = reshape(relative_diff, [], 1);

        rho_i     = sum(error_i_flat .^ 2);
        norm_z_sq = sum(relative_diff_flat .^ 2);

        N_i = neighbor_count(agent_idx);
        rho_bar_i = max((et_sigma * et_a * (1 - et_a * N_i) / N_i) * norm_z_sq, 1e-8);

        % 触发：||e_i||² > rho_bar_i
        if rho_i > rho_bar_i
            Zeta_last_trigger(:, agent_idx, :) = Zeta_vector(:, agent_idx, :);
            trigger_count(agent_idx) = 1;
        end
    end
end

%% Step 2: 从收敛后的 Zeta 提取聚合预测
Xi_all = P - Zeta_vector;   % [p_dim × AgentQuantity × M]

switch method
    case {'poe', 'gpoe', 'moe'}
        num1 = squeeze(Xi_all(1, :, :));
        den1 = squeeze(Xi_all(2, :, :));
        num2 = squeeze(Xi_all(3, :, :));
        den2 = squeeze(Xi_all(4, :, :));
        phi1 = num1 ./ max(den1, eps);
        phi2 = num2 ./ max(den2, eps);

    case 'bcm'
        num1 = squeeze(Xi_all(1, :, :));
        num2 = squeeze(Xi_all(2, :, :));
        den1 = squeeze(Xi_all(3, :, :));
        den2 = squeeze(Xi_all(4, :, :));
        prior_correction = (1 - AgentQuantity) / prior_var;
        den1_fused = den1 + prior_correction;
        den2_fused = den2 + prior_correction;
        phi1 = zeros(size(num1));
        phi2 = zeros(size(num2));
        mask1 = den1_fused > 1e-2;
        mask2 = den2_fused > 1e-2;
        phi1(mask1) = num1(mask1) ./ den1_fused(mask1);
        phi2(mask2) = num2(mask2) ./ den2_fused(mask2);

    case 'rbcm'
        num1  = squeeze(Xi_all(1, :, :));
        den1  = squeeze(Xi_all(2, :, :));
        beta1 = squeeze(Xi_all(3, :, :));
        num2  = squeeze(Xi_all(4, :, :));
        den2  = squeeze(Xi_all(5, :, :));
        beta2 = squeeze(Xi_all(6, :, :));
        den1_fused = den1 + (1 - beta1) / prior_var;
        den2_fused = den2 + (1 - beta2) / prior_var;
        phi1 = zeros(size(num1));
        phi2 = zeros(size(num2));
        mask1 = den1_fused > 1e-2;
        mask2 = den2_fused > 1e-2;
        phi1(mask1) = num1(mask1) ./ den1_fused(mask1);
        phi2(mask2) = num2(mask2) ./ den2_fused(mask2);
end

%% Step 3: 重建 MaskedGP
MaskedGP = cell(AgentQuantity, 1);
for AgentNr = 1:AgentQuantity
    Y_agent = [phi1(AgentNr, :); phi2(AgentNr, :)];
    MaskedGP{AgentNr} = LocalGP_MultiOutput(x_dim, y_dim, M, 1e-4, SigmaF, SigmaL);
    MaskedGP{AgentNr}.add_Alldata(InducingPoints_Coordinates, Y_agent);
end

end