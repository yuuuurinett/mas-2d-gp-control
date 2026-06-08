function [MaskedGP, Zeta_vector, Zeta_last_trigger, trigger_count] = gp_masked_aggregation_update( ...
    P, Zeta_vector, L, ...
    Kappa_P, AgentQuantity, NumInducingPoints, TimeStep, ...
    InducingPoints_Coordinates, SigmaF, SigmaL, x_dim, method, p_dim, ...
    Zeta_last_trigger, et_sigma, et_a, neighbor_count)

method    = lower(method);
M         = NumInducingPoints;
y_dim     = 2;
prior_var = SigmaF^2;

trigger_count = zeros(AgentQuantity, 1);

%% Step 1: DAC 积分（ode45，和原版一致，保证数值稳定）
if TimeStep > 0
    for InducingPointIdx = 1:M
        P_InducingPoint    = P(:, :, InducingPointIdx);
        Zeta_InducingPoint = Zeta_vector(:, :, InducingPointIdx);

        New_Consensus_Zeta_function = @(~, Zeta_ODE_Intern) Compute_New_Consensus_Derivative( ...
            Zeta_ODE_Intern, P_InducingPoint, L, Kappa_P, AgentQuantity, p_dim);

        [~, Zeta_ODE_Output] = ode45(New_Consensus_Zeta_function, ...
            [0, TimeStep], Zeta_InducingPoint(:));

        Zeta_vector(:, :, InducingPointIdx) = reshape(Zeta_ODE_Output(end,:)', ...
            p_dim, AgentQuantity);
    end

    %% Step 2: ET 触发判定（Dimarogonas 2012 公式11）
    % ET 控制的是"是否需要重建 MaskedGP"
    % 对每个 agent，把所有诱导点展平成一个大向量判断
    for agent_idx = 1:AgentQuantity
        % e_i = Zeta_last_trigger_i - Zeta_i
        error_i      = Zeta_last_trigger(:, agent_idx, :) - Zeta_vector(:, agent_idx, :);
        error_i_flat = reshape(error_i, [], 1);

        % z_i = sum_{j∈N_i}(Zeta_i - Zeta_j)
        neighbor_mask      = reshape(L(agent_idx,:) < 0, 1, AgentQuantity, 1);
        relative_diff      = sum((Zeta_vector(:,agent_idx,:) - Zeta_vector) .* neighbor_mask, 2);
        relative_diff_flat = reshape(relative_diff, [], 1);

        rho_i     = sum(error_i_flat .^ 2);
        norm_z_sq = sum(relative_diff_flat .^ 2);

        N_i       = neighbor_count(agent_idx);
        % rho_bar_i 加上下界（防止收敛后持续触发）和上界（防止阈值过大触发太少）
        rho_bar_i = (et_sigma * et_a * (1 - et_a * N_i) / N_i) * norm_z_sq;
        rho_bar_i = min(rho_bar_i, 1e-3);  % 上界
        rho_bar_i = max(rho_bar_i, 1e-8);  % 下界

        if rho_i > rho_bar_i
            Zeta_last_trigger(:, agent_idx, :) = Zeta_vector(:, agent_idx, :);
            trigger_count(agent_idx) = 1;
        end
    end
end

%% Step 3: 从 Zeta 提取聚合预测，重建 MaskedGP
Xi_all = P - Zeta_vector;

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

phi1(~isfinite(phi1)) = 0;
phi2(~isfinite(phi2)) = 0;

MaskedGP = cell(AgentQuantity, 1);
for AgentNr = 1:AgentQuantity
    Y_agent = [phi1(AgentNr, :); phi2(AgentNr, :)];
    MaskedGP{AgentNr} = LocalGP_MultiOutput(x_dim, y_dim, M, 1e-4, SigmaF, SigmaL);
    MaskedGP{AgentNr}.add_Alldata(InducingPoints_Coordinates, Y_agent);
end

end

%% 子函数：计算 DAC 动力学导数
function dZeta_dt = Compute_New_Consensus_Derivative(...
    Zeta_vec, P_Ref, L, Kappa, AgentQuantity, p_dim)

Zeta = reshape(Zeta_vec, p_dim, AgentQuantity);
dZeta_dt = Kappa * (P_Ref - Zeta) * L';
dZeta_dt = dZeta_dt(:);
end