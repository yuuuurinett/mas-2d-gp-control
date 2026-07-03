function dx_all = ET_MAS_GP_Leader_MultiAgent_DynamicFunction(...
	t,x_all,u_all,AgentQuantity,q_dim,SystemOrder)
dx_all = nan(size(x_all));
x_dim = q_dim * SystemOrder;
for AgentNr = 1:AgentQuantity
	x_i = x_all((AgentNr - 1) * x_dim + (1:x_dim));
	u_i = u_all((AgentNr - 1) * q_dim + (1:q_dim));
	dx_i = ET_MAS_GP_Leader_SingleAgent_DynamicFunction(t,x_i,u_i,q_dim,SystemOrder);
	dx_all((AgentNr - 1) * x_dim + (1:x_dim)) = dx_i;
end
end