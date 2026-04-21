function u_cell = Manipulator_2D_2DoF_get_u_cell(x_all_cell,phi_cell,f_hat_matrix, ...
	L1,L2,m1,m2)

AgentQuantity = numel(x_all_cell);
u_cell = cell(AgentQuantity,1);
for AgentNr = 1:AgentQuantity
	x_i = x_all_cell{AgentNr};
	phi_i = phi_cell{AgentNr};
	f_hat_i = f_hat_matrix(:,AgentNr);
	[h_i,g_i,~] = Manipulator_2D_2DoF_get_Dynamics_h_g_f(x_i,L1,L2,m1,m2);
	u_cell{AgentNr} = g_i \ (phi_i - h_i - f_hat_i);
end

end