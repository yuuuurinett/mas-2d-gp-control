function [Phi_Xi_vector, Zeta_vector] = gp_test_rbcm( ...
    AgentState_matrix, LocalGP_set, L, ...
    Kappa_P, AgentQuantity, Zeta_vector, TimeStep)

AgentQuantity = int32(AgentQuantity);

P_ReferenceSignal_matrix = zeros(6, AgentQuantity);
prior_variance = LocalGP_set{1}.SigmaF^2;

for AgentNr = 1:AgentQuantity
    [mu_n, var_n] = LocalGP_set{AgentNr}.predict(AgentState_matrix(:, AgentNr));
    beta_n  = 0.5 * (log(prior_variance) - log(var_n));
    beta_n1 = beta_n(1);
    beta_n2 = beta_n(2);
    var_n1  = var_n(1);
    var_n2  = var_n(2);

    P_ReferenceSignal_matrix(1, AgentNr) = AgentQuantity * beta_n1 * mu_n(1) / var_n1;
    P_ReferenceSignal_matrix(2, AgentNr) = AgentQuantity * beta_n1 / var_n1;
    P_ReferenceSignal_matrix(3, AgentNr) = AgentQuantity * beta_n1;
    P_ReferenceSignal_matrix(4, AgentNr) = AgentQuantity * beta_n2 * mu_n(2) / var_n2;
    P_ReferenceSignal_matrix(5, AgentNr) = AgentQuantity * beta_n2 / var_n2;
    P_ReferenceSignal_matrix(6, AgentNr) = AgentQuantity * beta_n2;
end

New_Consensus_Zeta_function = @(~, Zeta_vec) Compute_New_Consensus_Derivative( ...
    Zeta_vec, P_ReferenceSignal_matrix, L, Kappa_P, AgentQuantity);
[~, Zeta_Output] = ode45(New_Consensus_Zeta_function, [0, TimeStep], Zeta_vector(:));
Zeta_vector = reshape(Zeta_Output(end,:)', 6, AgentQuantity);

Xi = P_ReferenceSignal_matrix - Zeta_vector;

% dim1
Xi_num1  = Xi(1,:);   % β1·μ1/σ1²
Xi_den1  = Xi(2,:);   % β1/σ1²
Xi_beta1 = Xi(3,:);   % β1

% dim2  
Xi_num2  = Xi(4,:);   % β2·μ2/σ2²
Xi_den2  = Xi(5,:);   % β2/σ2²
Xi_beta2 = Xi(6,:);   % β2

denom1 = Xi_den1 + (1 - Xi_beta1) / prior_variance;
denom2 = Xi_den2 + (1 - Xi_beta2) / prior_variance;

% 完全 1:1 按照 rBCM 理论公式复刻
denom1 = Xi_den1 + (1 - Xi_beta1) / prior_variance;
denom2 = Xi_den2 + (1 - Xi_beta2) / prior_variance;

Phi_Xi_vector(1,:) = Xi_num1 ./ denom1;
Phi_Xi_vector(2,:) = Xi_num2 ./ denom2;



end

%%
function dZeta_dt = Compute_New_Consensus_Derivative( ...
    Zeta_vec, P_Ref, L, Kappa, AgentQuantity)

Zeta = reshape(Zeta_vec, 6, AgentQuantity);
dZeta_dt = Kappa * (P_Ref - Zeta) * L';
dZeta_dt = dZeta_dt(:);
end