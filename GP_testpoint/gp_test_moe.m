function [Phi_Xi_vector, Zeta_vector] = gp_test_moe( ...
    AgentState_matrix, LocalGP_set, L, ...
    Kappa_P, AgentQuantity, Zeta_vector, TimeStep)

P = zeros(4, AgentQuantity);  % 2 output dims x 2 moments each
for AgentNr = 1:AgentQuantity
    [mu_n, var_n] = LocalGP_set{AgentNr}.predict(AgentState_matrix(:, AgentNr));
    P(1, AgentNr) = AgentQuantity * mu_n(1);
    P(2, AgentNr) = AgentQuantity * mu_n(2);
    P(3, AgentNr) = AgentQuantity * (var_n(1) + mu_n(1)^2);
    P(4, AgentNr) = AgentQuantity * (var_n(2) + mu_n(2)^2);
end

New_Consensus_Zeta_function = @(~, Zeta_vec)  ...
    Compute_New_Consensus_Derivative( ...
    Zeta_vec, P, L, Kappa_P, AgentQuantity);

[~, Zeta_Output] = ode45(New_Consensus_Zeta_function, ...
    [0, TimeStep], Zeta_vector(:));
Zeta_vector = reshape(Zeta_Output(end,:)', 4, AgentQuantity);

Xi_matrix = P - Zeta_vector;

Phi_Xi_vector(1,:) = Xi_matrix(1,:) / AgentQuantity;
Phi_Xi_vector(2,:) = Xi_matrix(2,:) / AgentQuantity;
end

%% compute zeta_dot
function dZeta_dt = Compute_New_Consensus_Derivative(...
    Zeta_vec, P_Ref, L, Kappa, AgentQuantity)

Zeta = reshape(Zeta_vec, 4, AgentQuantity);
dZeta_dt = Kappa * (P_Ref - Zeta) * L';
dZeta_dt = dZeta_dt(:);

end