function [Phi_Xi_vector, Zeta_vector] = gp_rbcm( ...
    AgentState_matrix, LocalGP_set, ...
    AgentQuantity, Zeta_vector, TimeStep)

P_ReferenceSignal_matrix = zeros(6, AgentQuantity);
prior_variance= LocalGP_set{1}.SigmaF^2;
for AgentNr = 1:AgentQuantity
    [mu_n, var_n] = LocalGP_set{AgentNr}.predict(AgentState_matrix(:, AgentNr));
    %prior_variance= LocalGP_set{AgentNr}.SigmaF^2;
    beta_n = 0.5 * (log(prior_variance) - log(var_n));
    beta_n1=beta_n(1);
    beta_n2=beta_n(2);
    var_n1 = var_n(1);
    var_n2 = var_n(2);

    P_ReferenceSignal_matrix(1, AgentNr) = AgentQuantity * beta_n1 * mu_n(1) / var_n1;
    P_ReferenceSignal_matrix(2, AgentNr) = AgentQuantity * beta_n1 / var_n1;
    P_ReferenceSignal_matrix(3, AgentNr) = AgentQuantity * beta_n1;
    P_ReferenceSignal_matrix(4, AgentNr) = AgentQuantity * beta_n2 * mu_n(2) / var_n2;
    P_ReferenceSignal_matrix(5, AgentNr) = AgentQuantity * beta_n2 / var_n2;
    P_ReferenceSignal_matrix(6, AgentNr) = AgentQuantity * beta_n2;
end

New_Consensus_Zeta_function = @(~, Zeta_vec) Compute_New_Consensus_Derivative( ...
    Zeta_vec, P_ReferenceSignal_matrix, AgentQuantity);
[~, Zeta_Output] = ode45(New_Consensus_Zeta_function, [0, TimeStep], Zeta_vector(:));
Zeta_vector = reshape(Zeta_Output(end,:)', 6, AgentQuantity);

Xi_matrix = P_ReferenceSignal_matrix - Zeta_vector;
Phi_Xi_vector(1,:) = Xi_matrix(1,:) ./ (Xi_matrix(2,:) + (1 - Xi_matrix(3,:) / prior_variance));
Phi_Xi_vector(2,:) = Xi_matrix(4,:) ./ (Xi_matrix(5,:) + (1 - Xi_matrix(6,:) / prior_variance));
end

%%
function dZeta_dt = Compute_New_Consensus_Derivative(...
    Zeta_vec, P_Ref, L, Kappa, AgentQuantity)

Zeta = reshape(Zeta_vec, 4, AgentQuantity);
dZeta_dt = Kappa * (P_Ref - Zeta) * L';
dZeta_dt = dZeta_dt(:);

end