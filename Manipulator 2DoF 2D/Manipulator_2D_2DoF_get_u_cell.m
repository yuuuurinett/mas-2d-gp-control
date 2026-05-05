function u_cell = Manipulator_2D_2DoF_get_u_cell(x_all_cell,phi_cell,f_hat_matrix, ...
	L1,L2,m1,m2)

AgentQuantity = numel(x_all_cell);
u_cell = cell(AgentQuantity,1);
for AgentNr = 1:AgentQuantity
	x_i = x_all_cell{AgentNr};
	phi_i = phi_cell{AgentNr};
	f_hat_i = f_hat_matrix(:,AgentNr);

    if any(isnan(x_i)) || any(isinf(x_i)) || any(isnan(f_hat_i)) || any(isinf(f_hat_i))
        fprintf('🚨 致命错误拦截: Agent %d 在计算控制律前，状态或 GP 预测已变为 NaN/Inf！\n', AgentNr);
        % 强行输出 0 力矩，防止报错中断
        u_cell{AgentNr} = zeros(size(phi_i));
        continue; 
    end

	[h_i,g_i,~] = Manipulator_2D_2DoF_get_Dynamics_h_g_f(x_i,L1,L2,m1,m2);
    g_i_safe = g_i + 1e-6 * eye(size(g_i));

    try
        u_cell{AgentNr} = g_i_safe \ (phi_i - h_i - f_hat_i);
    catch
        fprintf('⚠️ 警告: Agent %d 的 g_i 矩阵发生奇异，已切换为伪逆 (pinv) 计算！\n', AgentNr);
        u_cell{AgentNr} = pinv(g_i_safe) * (phi_i - h_i - f_hat_i);
        u_cell{AgentNr} = max(-200, min(200, u_cell{AgentNr})); 
    end

	%u_cell{AgentNr} = g_i \ (phi_i - h_i - f_hat_i);
    %u_cell{AgentNr} = pinv(g_i) * (phi_i - h_i - f_hat_i);
end

end