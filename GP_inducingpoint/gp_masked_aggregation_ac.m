function MaskedGP = gp_masked_aggregation_ac( ...
    LocalGP_set, InducingPoints_Coordinates, SigmaF, SigmaL, ...
    x_dim, AgentQuantity, NumInducingPoints, method)
%% gp_masked_aggregation_ac

method    = lower(method);
M         = NumInducingPoints;
y_dim     = 2;
prior_var = SigmaF^2;

phi1 = zeros(AgentQuantity, M);  % AgentQuantity x M
phi2 = zeros(AgentQuantity, M);

for m = 1:M
    x_m = InducingPoints_Coordinates(:, m);

    % Collect local predictions from all agents
    mu_all  = zeros(AgentQuantity, y_dim);
    var_all = zeros(AgentQuantity, y_dim);
    for n = 1:AgentQuantity
        [mu_n, var_n]  = LocalGP_set{n}.predict(x_m);
        mu_all(n, :)  = mu_n';
        var_all(n, :) = var_n';
    end

    % Apply aggregation formula directly (no DAC needed)
for d = 1:y_dim    
    mu_d = mu_all(:, d);
    var_d = var_all(:, d);
        switch method
            case 'poe'
                prec = sum(1 ./ var_d);
                mu_fused = sum(mu_d ./ var_d) / prec;

            case 'gpoe'
                beta = max(eps, 0.5 * (log(prior_var) - log(var_d)));
                prec = sum(beta ./ var_d);
                mu_fused = sum(beta .* mu_d ./ var_d) / prec;

            case 'moe'
                mu_fused = mean(mu_d);

            case 'bcm'
                prec = sum(1 ./ var_d) + (1 - AgentQuantity) / prior_var;
                if prec > 1e-2
                    mu_fused = sum(mu_d ./ var_d) / prec;
                else
                    mu_fused = 0; 
                end
                %mu_fused = sum(mu_d ./ var_d) / prec;

            case 'rbcm'
                beta = max(eps, 0.5 * (log(prior_var) - log(var_d)));
                prec = sum(beta ./ var_d) + (1 - sum(beta)) / prior_var;
                if prec > 1e-2
                    mu_fused = sum(beta .* mu_d ./ var_d) / prec;
                else
                    mu_fused = 0; 
                end
        end

        if d == 1
            phi1(:, m) = mu_fused;
        else
            phi2(:, m) = mu_fused;
        end
end


%% Rebuild one MaskedGP per agent from fused inducing-point labels
MaskedGP = cell(AgentQuantity, 1);
for AgentNr = 1:AgentQuantity
    Y_agent = [phi1(AgentNr, :); phi2(AgentNr, :)];
    MaskedGP{AgentNr} = LocalGP_MultiOutput(x_dim, y_dim, M, 1e-4, SigmaF, SigmaL);
    MaskedGP{AgentNr}.add_Alldata(InducingPoints_Coordinates, Y_agent);
end
end