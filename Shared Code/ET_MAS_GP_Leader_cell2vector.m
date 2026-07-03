function a_all = ET_MAS_GP_Leader_cell2vector(a_cell, q_dim)
AgentQuantity = size(a_cell,1);
SystemOrder = size(a_cell,2);
x_dim = q_dim * SystemOrder;
a_all = nan(AgentQuantity * x_dim, 1);
for AgentNr = 1:AgentQuantity
	for SystemOrderNr = 1:SystemOrder
		q_ik = a_cell{AgentNr,SystemOrderNr};
		Index = (AgentNr - 1) * x_dim + (SystemOrderNr - 1) * q_dim + (1:q_dim);
		a_all(Index) = q_ik;
	end
end
end