function [mu_i_xi,var_i_xi,eta_i_xi] = Manipulator_2D_2DoF_LocalPrediction( ...
	x_i,AgentNr,LocalGP_set,beta,gamma,y_dim)
warning off;
[mu_i_xi, var_i_xi_set] = LocalGP_set{AgentNr}.predict(x_i);
warning on;
var_i_xi = var_i_xi_set(1);
eta_i_xi = sqrt(y_dim) * (sqrt(beta * var_i_xi) + gamma);

end