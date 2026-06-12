function [MaskedGP, Zeta_vector, Zeta_last_trigger, trigger_count] = gp_masked_aggregation_update( ...
    P, Zeta_vector, L, ...
    Kappa_P, AgentQuantity, NumInducingPoints, TimeStep, ...
    InducingPoints_Coordinates, SigmaF, SigmaL, x_dim, method, p_dim, ...
    Zeta_last_trigger, neighbor_count)

method    = lower(method);
M         = NumInducingPoints;
y_dim     = 2;
prior_var = SigmaF^2;

trigger_count = zeros(AgentQuantity, 1);

%% Step 1: DAC 积分（ode45，邻居项用Zeta_last_trigger，符合Kia 2014公式3）
if TimeStep > 0
    for InducingPointIdx = 1:M
        P_InducingPoint        = P(:, :, InducingPointIdx);
        Zeta_InducingPoint     = Zeta_vector(:, :, InducingPointIdx);
        Zeta_hat_InducingPoint = Zeta_last_trigger(:, :, InducingPointIdx);

        New_Consensus_Zeta_function = @(~, Zeta_ODE_Intern) Compute_New_Consensus_Derivative( ...
            Zeta_ODE_Intern, P_InducingPoint, Zeta_hat_InducingPoint, L, Kappa_P, AgentQuantity, p_dim);

        [~, Zeta_ODE_Output] = ode45(New_Consensus_Zeta_function, ...
            [0, TimeStep], Zeta_InducingPoint(:));

        Zeta_vector(:, :, InducingPointIdx) = reshape(Zeta_ODE_Output(end,:)', ...
            p_dim, AgentQuantity);
    end

    %% Step 2: ET触发判定（Kia 2014 CDC，公式17，针对无向连通图）
    % 触发条件：||e_tilde_i||² > (1/4d^i_out) * Σ_{j∈N_i}|x̂^i - x̂^j|² + (1/4d^i_out) * ε²
    % e_tilde_i = Zeta_last_trigger^i - Zeta^i（广播快照与当前值之差）
    % Pi 在 control task 中是时变的，因此使用 Kia 2014 动态 ET
    epsilon_i = 0.01;  % 兜底常数，防止 Zeno 行为
    for agent_i = 1:AgentQuantity
        d_i_out         = neighbor_count(agent_i);
        neighbor_mask_i = (L(agent_i,:) < 0);

        % e_tilde_i：广播快照与实时值之差（展平所有诱导点）
        e_tilde_i_flat = reshape(...
            Zeta_last_trigger(:,agent_i,:) - Zeta_vector(:,agent_i,:), [], 1);

        % 邻居广播值之差（Kia 2014公式17右侧第一项）
        neighbor_broadcast_diff_sq = sum(...
            (Zeta_last_trigger(:,agent_i,:) - Zeta_last_trigger(:,neighbor_mask_i,:)).^2, 'all');

        % 触发阈值（公式17）
        threshold_i = (1/(4*d_i_out)) * neighbor_broadcast_diff_sq + ...
                      (1/(4*d_i_out)) * epsilon_i^2;

        if sum(e_tilde_i_flat.^2) > threshold_i
            Zeta_last_trigger(:,agent_i,:) = Zeta_vector(:,agent_i,:);
            trigger_count(agent_i) = 1;
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

%% 子函数：计算 DAC 动力学导数（Kia 2014 公式3）
% ẋ^i = -α(x^i - u^i) - β Σ_j a_ij(x̂^i - x̂^j)
% 本地跟踪项用实时Zeta，邻居交互项用广播快照Zeta_last_trigger
function dZeta_dt = Compute_New_Consensus_Derivative(...
    Zeta_vec, P_Ref, Zeta_hat, L, Kappa, AgentQuantity, p_dim)

Zeta     = reshape(Zeta_vec, p_dim, AgentQuantity);
% 本地跟踪项：用实时Zeta
local_track = P_Ref - Zeta;
% 邻居交互项：用广播快照Zeta_hat（x̂^i - x̂^j）
neighbor_interact = Zeta_hat * L';
dZeta_dt = Kappa * (local_track - neighbor_interact);
dZeta_dt = dZeta_dt(:);
end