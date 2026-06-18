function P = recompute_P_single_agent(P, AgentNr, LocalGP_set, ...
    NumInducingPoints, InducingPoints_Coordinates, AgentQuantity, method)
% RECOMPUTE_P_SINGLE_AGENT
%   After agent AgentNr adds a new data point to its LocalGP, recompute
%   P(:, AgentNr, :) so that the DAC/AC state reflects the updated
%   posterior.  All other agents' columns are left untouched.
%
%   P      : p_dim x AgentQuantity x M  (modified in-place for AgentNr)
%   method : 'poe' | 'gpoe' | 'moe' | 'bcm' | 'rbcm'

method    = lower(method);
M         = NumInducingPoints;
prior_var = LocalGP_set{AgentNr}.SigmaF^2;

for m = 1:M
    x_m = InducingPoints_Coordinates(:, m);
    [mu_n, var_n] = LocalGP_set{AgentNr}.predict(x_m);

    mu1 = mu_n(1);  var1 = var_n(1);
    mu2 = mu_n(2);  var2 = var_n(2);

    switch method
        case 'poe'
            P(1, AgentNr, m) = AgentQuantity * mu1 / var1;
            P(2, AgentNr, m) = AgentQuantity / var1;
            P(3, AgentNr, m) = AgentQuantity * mu2 / var2;
            P(4, AgentNr, m) = AgentQuantity / var2;

        case 'gpoe'
            beta1 = max(eps, 0.5*(log(prior_var) - log(var1)));
            beta2 = max(eps, 0.5*(log(prior_var) - log(var2)));
            P(1, AgentNr, m) = AgentQuantity * beta1 * mu1 / var1;
            P(2, AgentNr, m) = AgentQuantity * beta1 / var1;
            P(3, AgentNr, m) = AgentQuantity * beta2 * mu2 / var2;
            P(4, AgentNr, m) = AgentQuantity * beta2 / var2;

        case 'moe'
            omega_n = 1.0 / AgentQuantity;
            P(1, AgentNr, m) = AgentQuantity * omega_n * mu1;
            P(2, AgentNr, m) = AgentQuantity * omega_n;
            P(3, AgentNr, m) = AgentQuantity * omega_n * mu2;
            P(4, AgentNr, m) = AgentQuantity * omega_n;

        case 'bcm'
            P(1, AgentNr, m) = AgentQuantity * mu1 / var1;
            P(2, AgentNr, m) = AgentQuantity * mu2 / var2;
            P(3, AgentNr, m) = AgentQuantity / var1;
            P(4, AgentNr, m) = AgentQuantity / var2;

        case 'rbcm'
            beta1 = max(eps, 0.5*(log(prior_var) - log(var1)));
            beta2 = max(eps, 0.5*(log(prior_var) - log(var2)));
            P(1, AgentNr, m) = AgentQuantity * beta1 * mu1 / var1;
            P(2, AgentNr, m) = AgentQuantity * beta1 / var1;
            P(3, AgentNr, m) = AgentQuantity * beta1;
            P(4, AgentNr, m) = AgentQuantity * beta2 * mu2 / var2;
            P(5, AgentNr, m) = AgentQuantity * beta2 / var2;
            P(6, AgentNr, m) = AgentQuantity * beta2;
    end
end
end