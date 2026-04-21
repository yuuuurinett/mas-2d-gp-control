function trigger_flag = Manipulator_2D_2DoF_DistributedET(AgentNr, ...
	r_matrix,e_cell,eta_underline_set,eta_aggregated_vector, ...
	MultiAgentSystem,Fl,xi,chi,vartheta_bar)
[q_dim,AgentQuantity] = size(r_matrix);
SystemOrder = size(e_cell,2);
x_dim = q_dim * SystemOrder;
%% Trigger
b_ii = MultiAgentSystem.Agent_Leader_Topology.ConnectionMatrix(AgentNr);
eta_aggregated_i = eta_aggregated_vector(AgentNr);
rho_i = (eta_aggregated_i + (1 - b_ii) * Fl)^2;
%% Trigger Threshold
z_i = nan(x_dim,1);
r_i = r_matrix(:,AgentNr);
z_i(1:q_dim) = r_i;
for SystemOrderNr = 1:SystemOrder-1
	z_i(SystemOrderNr * q_dim + (1:q_dim)) = e_cell{AgentNr, SystemOrderNr};
end
eta_underline_i = eta_underline_set(AgentNr);
rho_bar_i = xi^(-2) * max(z_i' * z_i - chi^(-2) * vartheta_bar ^ 2 / AgentQuantity, 0) + ...
	(eta_underline_i + (1 - b_ii) * Fl)^2;
%%
if rho_i > rho_bar_i
	trigger_flag = 1;
else
	trigger_flag = 0;
end

end