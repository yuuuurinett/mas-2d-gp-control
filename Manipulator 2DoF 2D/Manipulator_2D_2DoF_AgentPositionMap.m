function [plot_rx_set,plot_ry_set] = Manipulator_2D_2DoF_AgentPositionMap(...
	x_all, AgentQuantity)

x_all_matrix = reshape(x_all,[],AgentQuantity);
plot_rx_set = x_all_matrix(1,:)';
plot_ry_set = x_all_matrix(2,:)';

plot_rx_set = [plot_rx_set;plot_rx_set(1)];
plot_ry_set = [plot_ry_set;plot_ry_set(1)];
end