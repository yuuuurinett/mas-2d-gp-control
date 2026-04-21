function [phi_cell,r_matrix,e_cell] = Manipulator_2D_2DoF_ConsensusLaw(vartheta_cell,x_tilde_cell,x_l_r, ...
	MultiAgentSystem,c,lambda_set,s_r_cell)
%%
[AgentQuantity, SystemOrder] = size(vartheta_cell);
e_cell = cell(AgentQuantity, SystemOrder);
phi_cell = cell(AgentQuantity,1);
q_dim = numel(s_r_cell{1});
r_matrix = nan(q_dim,AgentQuantity);
%%
for AgentNr = 1:AgentQuantity
	b_ii = MultiAgentSystem.Agent_Leader_Topology.ConnectionMatrix(AgentNr);
	AgentNeighbourSet = MultiAgentSystem.Agent_Topology.NeighbourSet{AgentNr};
	AgentNeighbourQuantity = numel(AgentNeighbourSet);
	r_i = 0;
	for SystemOrderNr = 1:SystemOrder
		e_ik = b_ii * vartheta_cell{AgentNr, SystemOrderNr};
		x_tilde_ik = x_tilde_cell{AgentNr, SystemOrderNr};
		for NeighborNr = 1:AgentNeighbourQuantity
			NeighborAgentNr = AgentNeighbourSet(NeighborNr);
			a_ij = MultiAgentSystem.Agent_Topology.AdjacencyMatrix(AgentNr, NeighborAgentNr);
			x_tilde_jk = x_tilde_cell{NeighborAgentNr, SystemOrderNr};
			e_ik = e_ik + a_ij * (x_tilde_ik - x_tilde_jk);
		end
		e_cell{AgentNr, SystemOrderNr} = e_ik;
		lambda_k = lambda_set(SystemOrderNr);
		r_i = r_i + lambda_k * e_ik;
	end
	r_matrix(:,AgentNr) = r_i;
	s_r_i = s_r_cell{AgentNr};
	phi_cell{AgentNr} = -c * r_i + s_r_i + b_ii * x_l_r;
end

end