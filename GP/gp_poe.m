function [Phi_Xi_vector, Zeta_vector, Xi_diff] = gp_poe( ...
    AgentState_matrix, LocalGP_set, L, ...
    Kappa_P, AgentQuantity, Zeta_vector, TimeStep)


P_ReferenceSignal = zeros(4, AgentQuantity);
for AgentNr = 1:AgentQuantity
    [mu_n, var_n] = LocalGP_set{AgentNr}.predict(AgentState_matrix(:, AgentNr));
    var_n1 = var_n(1);
    var_n2 = var_n(2);
    P_ReferenceSignal(1, AgentNr) = AgentQuantity * mu_n(1) / var_n1;
    P_ReferenceSignal(2, AgentNr) = AgentQuantity * mu_n(2) / var_n2;
    P_ReferenceSignal(3, AgentNr) = AgentQuantity / var_n1;
    P_ReferenceSignal(4, AgentNr) = AgentQuantity / var_n2;
end

New_Consensus_Zeta_function = @(~, Zeta_vec)  ...
    Compute_New_Consensus_Derivative( ...
    Zeta_vec, P_ReferenceSignal, L, Kappa_P, AgentQuantity);

[~, Zeta_Output] = ode45(New_Consensus_Zeta_function, ...
    [0, TimeStep], Zeta_vector(:));
Zeta_vector = reshape(Zeta_Output(end,:)', 4, AgentQuantity);


Xi_matrix = P_ReferenceSignal - Zeta_vector;

Phi_Xi_vector(1,:) = Xi_matrix(1,:) ./ Xi_matrix(3,:);
Phi_Xi_vector(2,:) = Xi_matrix(2,:) ./ Xi_matrix(4,:);

Xi_diff = sum(pdist(Xi_matrix', 'euclidean'));
end

%%
function dZeta_dt = Compute_New_Consensus_Derivative(...
    Zeta_vec, P_Ref, L, Kappa, AgentQuantity)

Zeta = reshape(Zeta_vec, 4, AgentQuantity);
dZeta_dt = Kappa * (P_Ref - Zeta) * L';
dZeta_dt = dZeta_dt(:);

end