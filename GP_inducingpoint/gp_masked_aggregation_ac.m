function[MaskedGP, ac_consensus_trigger_count] = gp_masked_aggregation_ac( ...
    LocalGP_set, InducingPoints_Coordinates, SigmaF, SigmaL, ...
    x_dim, AgentQuantity, NumInducingPoints, method)
%% gp_masked_aggregation_ac
% IP-AC：纯共识迭代 + Dimarogonas 2012 ET（公式11）
% 动力学：dXi/dt = -L * Xi（用广播快照 Xi_last_trigger）
% 触发条件：||e_i||² > σ_i * a * (1-a*N_i)/N_i * ||z_i||²
% z_i 用实时 Xi（Dimarogonas 2012 论文）

method    = lower(method);
M         = NumInducingPoints;
y_dim     = 2;
prior_var = SigmaF^2;
p_dim     = 2 * y_dim;

%% 1. 拓扑参数
% 用循环无向图（和dataset实验一致）
% 这里直接用传入的LocalGP_set推断AgentQuantity
MultiAgentSystem = Manipulator_2D_2DoF_SetMASTopology(AgentQuantity, 1);
L = MultiAgentSystem.Agent_Topology.LaplacianMatrix;
N_degree = sum(L < 0, 2);
N_max    = max(N_degree);
sigma_i  = 0.5;    % Dimarogonas 2012 设计参数
a_param  = 0.5 / N_max;

%% 2. 预计算诱导点上的局部预测，构建Pi矩阵
Pi = zeros(p_dim, AgentQuantity, M);
for n = 1:AgentQuantity
    for m = 1:M
        x_m = InducingPoints_Coordinates(:, m);
        [mu_n, var_n] = LocalGP_set{n}.predict(x_m);
        for d = 1:y_dim
            vs = max(var_n(d), 1e-6);
            b  = max(min(0.5*(log(prior_var)-log(vs)), 10), eps);
            switch method
                case 'poe'
                    Pi(2*d-1,n,m) = AgentQuantity * mu_n(d) / vs;
                    Pi(2*d,  n,m) = AgentQuantity / vs;
                case 'gpoe'
                    Pi(2*d-1,n,m) = AgentQuantity * b * mu_n(d) / vs;
                    Pi(2*d,  n,m) = AgentQuantity * b / vs;
                case 'moe'
                    Pi(2*d-1,n,m) = AgentQuantity * mu_n(d);
                    Pi(2*d,  n,m) = AgentQuantity * (vs + mu_n(d)^2);
                case 'bcm'
                    Pi(2*d-1,n,m) = AgentQuantity * mu_n(d) / vs;
                    Pi(2*d,  n,m) = AgentQuantity / vs - (AgentQuantity-1)/prior_var;
                case 'rbcm'
                    Pi(2*d-1,n,m) = AgentQuantity * b * mu_n(d) / vs;
                    Pi(2*d,  n,m) = AgentQuantity * b / vs + (1 - AgentQuantity*b)/prior_var;
            end
        end
    end
end

%% 3. AC迭代共识 + Dimarogonas 2012 ET
Kappa_P = 10; t_step = 0.01; max_iter = 3000;
Xi              = Pi;   % 实时状态 x_i(t)，初始值 = Pi
Xi_last_trigger = Pi;   % 广播快照 x̂_i(t)，初始值 = Pi
ac_consensus_trigger_count = zeros(AgentQuantity, M);
for iter = 1:max_iter
    Xi_prev = Xi;

    % 动力学向量化：xdot = -Kappa_P * L * x̂（用广播快照，零阶保持）
    Xi_hat_perm = permute(Xi_last_trigger, [2 1 3]);
    Xi_hat_2d   = reshape(Xi_hat_perm, AgentQuantity, p_dim*M);
    L_Xi_hat_2d = L * Xi_hat_2d;
    L_Xi_hat    = permute(reshape(L_Xi_hat_2d, AgentQuantity, p_dim, M), [2 1 3]);
    Xi = Xi - t_step * Kappa_P * L_Xi_hat;

    % ET触发条件（Dimarogonas 2012，公式11，point-wise向量化）
    % e_i = x̂_i - x_i（广播快照与实时值之差）
    % z_i = Σ_{j∈N_i}(x_i - x_j)（实时状态的邻居差）
    for agent_i = 1:AgentQuantity
        neighbor_mask_i = (L(agent_i,:) < 0);
        N_i     = N_degree(agent_i);
        coeff_i = sigma_i * a_param * (1 - a_param*N_i) / N_i;

        % e_tilde_all：[p_dim x 1 x M]
        e_tilde_all = Xi_last_trigger(:,agent_i,:) - Xi(:,agent_i,:);
        % z_all：实时邻居差，[p_dim x 1 x M]
        z_all = sum(Xi(:,agent_i,:) - Xi(:,neighbor_mask_i,:), 2);

        e_sq_all      = squeeze(sum(e_tilde_all.^2, 1));  % [M x 1]
        z_sq_all      = squeeze(sum(z_all.^2,       1));  % [M x 1]
        threshold_all = coeff_i * z_sq_all;               % [M x 1]

        trigger_mask = (e_sq_all > threshold_all) | (iter == 1);
        Xi_last_trigger(:,agent_i,trigger_mask) = Xi(:,agent_i,trigger_mask);
        ac_consensus_trigger_count(agent_i, trigger_mask) = ...
        ac_consensus_trigger_count(agent_i, trigger_mask) + 1;
    end

    if max(abs(Xi(:) - Xi_prev(:))) < 1e-5, break; end
end
% fprintf('[IP-AC control] 收敛步数:%d\n', iter);  % suppress batch output

%% 4. 提取phi并重建MaskedGP
MaskedGP = cell(AgentQuantity, 1);
for agent_i = 1:AgentQuantity
    phi = zeros(y_dim, M);
    for d = 1:y_dim
        xi1 = squeeze(Xi(2*d-1, agent_i, :))';
        xi2 = squeeze(Xi(2*d,   agent_i, :))';
        if ismember(method, {'gpoe','poe','bcm','rbcm'})
            phi(d,:) = xi1 ./ max(abs(xi2), 1e-4);
        else
            phi(d,:) = xi1 / AgentQuantity;
        end
    end
    MaskedGP{agent_i} = LocalGP_MultiOutput(x_dim, y_dim, M, 1e-4, SigmaF, SigmaL);
    MaskedGP{agent_i}.add_Alldata(InducingPoints_Coordinates, phi);
end
end