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

for m = 1:M
    x_m = InducingPoints_Coordinates(:, m);
    [mu_n, var_n] = LocalGP_set{AgentNr}.predict(x_m);

    for d = 1:y_dim
        mu = mu_n(d);
        vs = max(var_n(d), 1e-8);
        beta = 0.5 * (log(prior_var) - log(vs));
        beta = max(min(beta, 10), eps);

        switch method
            case 'poe'
                P(2*d-1, AgentNr, m) = AgentQuantity * mu / vs;
                P(2*d,   AgentNr, m) = AgentQuantity / vs;

            case 'gpoe'
                P(2*d-1, AgentNr, m) = AgentQuantity * beta * mu / vs;
                P(2*d,   AgentNr, m) = AgentQuantity * beta / vs;

            case 'moe'
                P(2*d-1, AgentNr, m) = AgentQuantity * mu;
                P(2*d,   AgentNr, m) = AgentQuantity * (vs + mu^2);

            case 'bcm'
                P(2*d-1, AgentNr, m) = AgentQuantity * mu / vs;
                P(2*d,   AgentNr, m) = AgentQuantity / vs - ...
                    (AgentQuantity - 1) / prior_var;

            case 'rbcm'
                P(2*d-1, AgentNr, m) = AgentQuantity * beta * mu / vs;
                P(2*d,   AgentNr, m) = AgentQuantity * beta / vs + ...
                    (1 - AgentQuantity * beta) / prior_var;

            otherwise
                error('Unknown aggregation method: %s', method);
        end
    end
end
end
