function P = recompute_P_single_agent(P, AgentNr, LocalGP_set, ...
    NumInducingPoints, InducingPoints_Coordinates, AgentQuantity, method)
% RECOMPUTE_P_SINGLE_AGENT
% Recompute only P(:,AgentNr,:) after an online LocalGP update.
% Encoding is consistent with gp_masked_aggregation_ac.m:
%   p_dim = 4 for all methods: [num1; den_or_aux1; num2; den_or_aux2].

method = lower(method);
M = NumInducingPoints;
prior_var = LocalGP_set{AgentNr}.SigmaF^2;
y_dim = 2;
aggregation_cfg = control_aggregation_parameters();
posterior_var_floor = aggregation_cfg.posterior_var_floor;
rbcm_beta_max = aggregation_cfg.rbcm_beta_max;
bcm_prior_scale = aggregation_cfg.bcm_prior_scale;

for m = 1:M
    x_m = InducingPoints_Coordinates(:, m);
    [mu_n, var_n] = LocalGP_set{AgentNr}.predict(x_m);

    for d = 1:y_dim
        % 均值截断：防止online learning初期GP外推产生离谱预测值
        mu = max(-30, min(30, mu_n(d)));
        % 方差下限：1e-8对GP而言太小，容易在数据点重合时让1/vs爆炸
        % 提高到1e-3，符合物理上GP不确定性的合理下限
        vs = max(var_n(d), posterior_var_floor);
        raw_beta = 0.5 * (log(prior_var) - log(vs));
        beta_gpoe = max(min(raw_beta,aggregation_cfg.gpoe_beta_max),eps);
        beta_rbcm = min(max(raw_beta, eps), rbcm_beta_max);

        switch method
            case 'poe'
                P(2*d-1, AgentNr, m) = AgentQuantity * mu / vs;
                P(2*d,   AgentNr, m) = AgentQuantity / vs;

            case 'gpoe'
                P(2*d-1, AgentNr, m) = AgentQuantity * beta_gpoe * mu / vs;
                P(2*d,   AgentNr, m) = AgentQuantity * beta_gpoe / vs;

            case 'moe'
                P(2*d-1, AgentNr, m) = AgentQuantity * mu;
                P(2*d,   AgentNr, m) = AgentQuantity * (vs + mu^2);

            case 'bcm'
                P(2*d-1, AgentNr, m) = AgentQuantity * mu / vs;
                P(2*d,   AgentNr, m) = AgentQuantity / vs - ...
                    bcm_prior_scale * (AgentQuantity - 1) / prior_var;

            case 'rbcm'
                P(2*d-1, AgentNr, m) = AgentQuantity * beta_rbcm * mu / vs;
                P(2*d,   AgentNr, m) = AgentQuantity * beta_rbcm / vs + ...
                    (1 - AgentQuantity * beta_rbcm) / prior_var;

            otherwise
                error('Unknown aggregation method: %s', method);
        end

        % 最后一道防线：限制P值幅度，防止任何残余的数值问题传播到consensus层
        P(2*d-1, AgentNr, m) = max(-1e4, min(1e4, P(2*d-1, AgentNr, m)));
        P(2*d,   AgentNr, m) = max(-1e4, min(1e4, P(2*d,   AgentNr, m)));
    end
end
end
