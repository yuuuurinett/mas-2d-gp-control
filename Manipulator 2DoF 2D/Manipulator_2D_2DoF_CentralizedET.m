function trigger_flag = Manipulator_2D_2DoF_CentralizedET( ...
	r_matrix,e_cell,eta_aggregated_vector, ...
	MultiAgentSystem,Fl,xi,chi,vartheta_bar)
[q_dim,AgentQuantity] = size(r_matrix);
SystemOrder = size(e_cell,2);
%%
B = diag(MultiAgentSystem.Agent_Leader_Topology.ConnectionMatrix);
rho = norm((eye(AgentQuantity) - B) * ones(AgentQuantity,1) * Fl + eta_aggregated_vector);
%%
epsilon = nan((SystemOrder - 1) * q_dim * AgentQuantity,1);
for SystemOrderNr = 1:SystemOrder - 1
	epsilon_k = nan(q_dim * AgentQuantity,1);
	for AgentNr = 1:AgentQuantity
		epsilon_k((AgentNr - 1) * q_dim + (1:q_dim)) = e_cell{AgentNr,SystemOrderNr};
	end
	epsilon((SystemOrderNr - 1) * AgentQuantity * q_dim + (1:(AgentQuantity * q_dim))) = epsilon_k;
end
r = nan(AgentQuantity * q_dim,1);
for AgentNr = 1:AgentQuantity
	r_i = r_matrix(:,AgentNr);
	r((AgentNr - 1) * q_dim + (1:q_dim)) = r_i;
end
z = [r;epsilon];
rho_bar = xi^(-1) * max(norm(z), chi^(-1) * vartheta_bar);
%%
if rho > rho_bar
	trigger_flag = 1;
else
	trigger_flag = 0;
end

end