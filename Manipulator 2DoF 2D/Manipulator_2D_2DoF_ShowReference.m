function Manipulator_2D_2DoF_ShowReference( ...
	t_set,xl_set,s_all_set,AgentQuantity,L1,L2)

x_dim = size(xl_set,1);

FigureObj_Trajectory = figure('Name','Reference');
AxesObj_Trajectory = axes(FigureObj_Trajectory);
axis(AxesObj_Trajectory,[-(L1+L2),(L1+L2),-(L1+L2),(L1+L2)]);
hold(AxesObj_Trajectory,'on');
plot(AxesObj_Trajectory,xl_set(1,:),xl_set(2,:),'k-.');
for AgentNr = 1:AgentQuantity
	plot(AxesObj_Trajectory, ...
		xl_set(1,:) + s_all_set((AgentNr - 1) * x_dim + 1,:), ...
		xl_set(2,:) + s_all_set((AgentNr - 1) * x_dim + 2,:), ...
		'--');
end
[plot_rx_set,plot_ry_set] = Manipulator_2D_2DoF_AgentPositionMap(...
	kron(ones(AgentQuantity,1),xl_set(:,1)) + s_all_set(:,1), AgentQuantity);
plot(AxesObj_Trajectory,plot_rx_set,plot_ry_set,'r.-');
for t_Nr = 2:numel(t_set)
	[plot_rx_set,plot_ry_set] = Manipulator_2D_2DoF_AgentPositionMap(...
		kron(ones(AgentQuantity,1),xl_set(:,t_Nr)) + s_all_set(:,t_Nr), AgentQuantity);
	AxesObj_Trajectory.Children(1).XData = plot_rx_set;
	AxesObj_Trajectory.Children(1).YData = plot_ry_set;
	drawnow;
	pause(0.001);
end

end