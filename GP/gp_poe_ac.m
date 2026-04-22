function [Phi_Xi_vector] = gp_poe_ac( ...
    AgentState_matrix, LocalGP_set, AgentQuantity)

P = zeros(4, AgentQuantity);
for AgentNr = 1:AgentQuantity
    [mu_n, var_n] = LocalGP_set{AgentNr}.predict(AgentState_matrix(:, AgentNr));
    P(1, AgentNr) = AgentQuantity * mu_n(1) / var_n(1);
    P(2, AgentNr) = AgentQuantity * mu_n(2) / var_n(2);
    P(3, AgentNr) = AgentQuantity / var_n(1);
    P(4, AgentNr) = AgentQuantity / var_n(2);
end

Xi = repmat(mean(P, 2), 1, AgentQuantity); 

Phi_Xi_vector(1,:) = Xi(1,:) ./ Xi(3,:);
Phi_Xi_vector(2,:) = Xi(2,:) ./ Xi(4,:);
end