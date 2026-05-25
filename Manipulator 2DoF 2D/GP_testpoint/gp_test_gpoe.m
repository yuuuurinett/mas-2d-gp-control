function [Phi_Xi_vector, Zeta_vector] = gp_test_gpoe( ...
    AgentState_matrix, LocalGP_set, L, ...
    Kappa_P, AgentQuantity, Zeta_vector, TimeStep)

prior_var = LocalGP_set{1}.SigmaF^2;
P = zeros(4, AgentQuantity);
for AgentNr = 1:AgentQuantity
    [mu_n, var_n] = LocalGP_set{AgentNr}.predict(AgentState_matrix(:, AgentNr));
    beta1 = max(eps, 0.5 * (log(prior_var) - log(var_n(1))));
    beta2 = max(eps, 0.5 * (log(prior_var) - log(var_n(2))));
    P(1, AgentNr) = AgentQuantity * beta1 * mu_n(1) / var_n(1);
    P(2, AgentNr) = AgentQuantity * beta2 * mu_n(2) / var_n(2);
    P(3, AgentNr) = AgentQuantity * beta1 / var_n(1);
    P(4, AgentNr) = AgentQuantity * beta2 / var_n(2);
end

New_Consensus_Zeta_function = @(~, Zeta_vec)  ...
    Compute_New_Consensus_Derivative( ...
    Zeta_vec, P, L, Kappa_P, AgentQuantity);

[~, Zeta_Output] = ode45(New_Consensus_Zeta_function, ...
    [0, TimeStep], Zeta_vector(:));
Zeta_vector = reshape(Zeta_Output(end,:)', 4, AgentQuantity);

Xi = P - Zeta_vector;
Phi_Xi_vector(1,:) = Xi(1,:) ./ Xi(3,:);
Phi_Xi_vector(2,:) = Xi(2,:) ./ Xi(4,:);
end

%% compute zeta_dot
function dZeta_dt = Compute_New_Consensus_Derivative(...
    Zeta_vec, P_Ref, L, Kappa, AgentQuantity)

Zeta = reshape(Zeta_vec, 4, AgentQuantity);
dZeta_dt = Kappa * (P_Ref - Zeta) * L';
dZeta_dt = dZeta_dt(:);

end