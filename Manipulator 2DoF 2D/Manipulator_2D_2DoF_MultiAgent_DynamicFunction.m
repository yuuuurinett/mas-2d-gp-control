function dx_all = Manipulator_2D_2DoF_MultiAgent_DynamicFunction(...
	t,x_all,u_cell,L1,L2,m1,m2)
dx_all = nan(size(x_all));
AgentQuantity = numel(u_cell);
q_dim = numel(u_cell{1});
SystemOrder = numel(x_all) / AgentQuantity / q_dim;
x_dim = q_dim * SystemOrder;
for AgentNr = 1:AgentQuantity
	x_i = x_all((AgentNr - 1) * x_dim + (1:x_dim));
	u_i = u_cell{AgentNr};
	dx_i = Manipulator_2D_2DoF_DynamicFunc(t,x_i,u_i,L1,L2,m1,m2);
	dx_all((AgentNr - 1) * x_dim + (1:x_dim)) = dx_i;
end
end