function dx_all = multi_agent_dynamics(t, x_all, ...
    Phi_Xi_vector, AdjMatrix, NeighbourSet, B, lambda_set, Kappa_C, ...
    AgentQuantity, SystemOrder, q_dim, L1, L2, m1, m2)

x_dim = q_dim * SystemOrder;
AgentState_matrix = reshape(x_all, x_dim, AgentQuantity);

LeaderState_vector = [sin(0.5*t); cos(0.5*t); 0.5*cos(0.5*t); -0.5*sin(0.5*t)];

% leader 加速度前馈 r̈_l
x_l_r = [-0.25*sin(0.5*t); -0.25*cos(0.5*t)];

dx_all = zeros(size(x_all));

for AgentNr = 1:AgentQuantity
    [Epsilon_vector, r_n] = consensus_error( ...
        AgentNr, AgentState_matrix, LeaderState_vector, ...
        AdjMatrix, NeighbourSet, B, lambda_set, SystemOrder, q_dim);

    x_i = AgentState_matrix(:, AgentNr);

    [h_i, g_i, ~] = Manipulator_2D_2DoF_get_Dynamics_h_g_f(x_i, L1, L2, m1, m2);

    % 加入 leader 加速度前馈
    Nu_i = controller(r_n, Epsilon_vector, Kappa_C, lambda_set, SystemOrder) + x_l_r;

    u_i = g_i \ (Nu_i - h_i - Phi_Xi_vector(:, AgentNr));

    dx_i = Manipulator_2D_2DoF_DynamicFunc(t, x_i, u_i, L1, L2, m1, m2);
    dx_all((AgentNr-1)*x_dim + (1:x_dim)) = dx_i;
end
end