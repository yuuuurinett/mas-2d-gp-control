function [mu_cell,var_matrix,eta_matrix] = ET_MAS_GP_Leader_LocalPrediction(AgentNr, ...
	x_matrix,LocalGP_set,beta,gamma,y_dim,mu_cell,var_matrix,eta_matrix)

x_i = x_matrix(:,AgentNr);
[mu_i_xi, var_i_xi] = LocalGP_set{AgentNr}.predict(x_i);
eta_i_xi = sqrt(y_dim) * (sqrt(beta * var_i_xi(1)) + gamma);
mu_cell{AgentNr, AgentNr} = mu_i_xi;
var_matrix(AgentNr, AgentNr) = var_i_xi(1);
eta_matrix(AgentNr, AgentNr) = eta_i_xi;
end