function [P, p_dim] = gp_masked_aggregation_init( ...
    LocalGP_set, AgentQuantity, ...
    NumInducingPoints, InducingPoints_Coordinates, method)
%    InducingPoints_Coordinates : x_dim x M
%    P      : p_dim x AgentQuantity x M
%    p_dim  : dimension of P's first axis (4 or 6)

method = lower(method);
M = NumInducingPoints;
prior_var = LocalGP_set{1}.SigmaF^2;  % sigma_0^²

switch method
    case {'poe', 'gpoe', 'moe', 'bcm', 'rbcm'}
        % Unified layout: [numerator_1; denominator_1; ...].  Scaling each
        % local contribution by AgentQuantity makes the DAC average equal
        % to the centralized sum required by the product-based methods.
        p_dim = 4;
    otherwise
        error('Unknown aggregation method: %s. Choose from poe/gpoe/moe/bcm/rbcm.', method);
end

P = zeros(p_dim, AgentQuantity, M);

for AgentNr = 1:AgentQuantity
    for InducingPointIdx = 1:M
        x_m = InducingPoints_Coordinates(:, InducingPointIdx);  % x_dim x 1
        [mu_n, var_n] = LocalGP_set{AgentNr}.predict(x_m);  % 2x1, 2x1

        mu1  = max(-30, min(30, mu_n(1)));
        var1 = max(var_n(1), 1e-3);
        mu2  = max(-30, min(30, mu_n(2)));
        var2 = max(var_n(2), 1e-3);

        switch method

            case 'poe'                        
                P(1, AgentNr, InducingPointIdx) = AgentQuantity * mu1 / var1;
                P(2, AgentNr, InducingPointIdx) = AgentQuantity / var1;
                P(3, AgentNr, InducingPointIdx) = AgentQuantity * mu2 / var2;
                P(4, AgentNr, InducingPointIdx) = AgentQuantity / var2;

            case 'gpoe'        
                beta1 = max(eps, 0.5 * (log(prior_var) - log(var1)));
                beta2 = max(eps, 0.5 * (log(prior_var) - log(var2)));
                P(1, AgentNr, InducingPointIdx) = AgentQuantity * beta1 * mu1 / var1;
                P(2, AgentNr, InducingPointIdx) = AgentQuantity * beta1 / var1;
                P(3, AgentNr, InducingPointIdx) = AgentQuantity * beta2 * mu2 / var2;
                P(4, AgentNr, InducingPointIdx) = AgentQuantity * beta2 / var2;

            case 'moe'
                P(1, AgentNr, InducingPointIdx) = AgentQuantity * mu1;
                P(2, AgentNr, InducingPointIdx) = AgentQuantity * (var1 + mu1^2);
                P(3, AgentNr, InducingPointIdx) = AgentQuantity * mu2;
                P(4, AgentNr, InducingPointIdx) = AgentQuantity * (var2 + mu2^2);

            case 'bcm'
                P(1, AgentNr, InducingPointIdx) = AgentQuantity * mu1 / var1;
                P(2, AgentNr, InducingPointIdx) = AgentQuantity / var1 - ...
                    (AgentQuantity - 1) / prior_var;
                P(3, AgentNr, InducingPointIdx) = AgentQuantity * mu2 / var2;
                P(4, AgentNr, InducingPointIdx) = AgentQuantity / var2 - ...
                    (AgentQuantity - 1) / prior_var;

            case 'rbcm'
                beta1 = max(eps, 0.5 * (log(prior_var) - log(var1)));
                beta2 = max(eps, 0.5 * (log(prior_var) - log(var2)));
                P(1, AgentNr, InducingPointIdx) = AgentQuantity * beta1 * mu1 / var1;
                P(2, AgentNr, InducingPointIdx) = AgentQuantity * beta1 / var1 + ...
                    (1 - AgentQuantity * beta1) / prior_var;
                P(3, AgentNr, InducingPointIdx) = AgentQuantity * beta2 * mu2 / var2;
                P(4, AgentNr, InducingPointIdx) = AgentQuantity * beta2 / var2 + ...
                    (1 - AgentQuantity * beta2) / prior_var;
        end
    end
end
end
