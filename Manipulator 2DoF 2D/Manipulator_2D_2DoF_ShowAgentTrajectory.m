function Manipulator_2D_2DoF_ShowAgentTrajectory( ...
	t_set,xl_set,s_all_set,x_all_set, ...
	L1,L2)
x_dim = size(xl_set,1);
AgentQuantity = size(s_all_set,1) / x_dim;
%%
FigureObj_Trajectory = figure('Name','Trajectory');
AxesObj_Trajectory = axes(FigureObj_Trajectory);
axis(AxesObj_Trajectory,[-(L1+L2),(L1+L2),-(L1+L2),(L1+L2)]);
hold(AxesObj_Trajectory,'on');
plot(AxesObj_Trajectory,xl_set(1,:),xl_set(2,:),'k-.');
%%
Agent_ReferenceTrajectoryObj_cell = cell(AgentQuantity,1);
Agent_StateTrajectoryObj_cell = cell(AgentQuantity,1);
% t_Nr = 1
x_all = x_all_set(:,1);
x_all_cell = ET_MAS_GP_Leader_vector2cell(x_all, AgentQuantity, 1);
for AgentNr = 1:AgentQuantity
	Agent_ReferenceTrajectoryObj_cell{AgentNr} = plot(AxesObj_Trajectory, ...
		xl_set(1,1) + s_all_set((AgentNr - 1) * x_dim + 1,1), ...
		xl_set(2,1) + s_all_set((AgentNr - 1) * x_dim + 2,1), ...
		'--');

	x_i = x_all_cell{AgentNr};
	Agent_StateTrajectoryObj_cell{AgentNr} = plot(AxesObj_Trajectory, ...
		x_i(1),x_i(2),'-');
end
[plot_rx_set,plot_ry_set] = Manipulator_2D_2DoF_AgentPositionMap(...
	x_all, AgentQuantity);
Agent_ActualFormationObj = plot(AxesObj_Trajectory, ...
	plot_rx_set,plot_ry_set,'r.-');
[plot_rx_ref_set,plot_ry_ref_set] = Manipulator_2D_2DoF_AgentPositionMap(...
	kron(ones(AgentQuantity,1),xl_set(:,1)) + s_all_set(:,1), AgentQuantity);
Agent_ReferenceFormationObj = plot(AxesObj_Trajectory, ...
	plot_rx_ref_set,plot_ry_ref_set,'b.--');
% t_Nr
for t_Nr = 2:numel(t_set)
	x_all = x_all_set(:,t_Nr);
	x_all_cell = ET_MAS_GP_Leader_vector2cell(x_all, AgentQuantity, 1);
	for AgentNr = 1:AgentQuantity
		Agent_ReferenceTrajectoryObj_cell{AgentNr}.XData = [ ...
			Agent_ReferenceTrajectoryObj_cell{AgentNr}.XData, ...
			xl_set(1,t_Nr) + s_all_set((AgentNr - 1) * x_dim + 1,t_Nr)];
		Agent_ReferenceTrajectoryObj_cell{AgentNr}.YData = [ ...
			Agent_ReferenceTrajectoryObj_cell{AgentNr}.YData, ...
			xl_set(2,t_Nr) + s_all_set((AgentNr - 1) * x_dim + 2,t_Nr)];

		x_i = x_all_cell{AgentNr};
		Agent_StateTrajectoryObj_cell{AgentNr}.XData = [ ...
			Agent_StateTrajectoryObj_cell{AgentNr}.XData,x_i(1)];
		Agent_StateTrajectoryObj_cell{AgentNr}.YData = [ ...
			Agent_StateTrajectoryObj_cell{AgentNr}.YData,x_i(2)];
	end
	[plot_rx_set,plot_ry_set] = Manipulator_2D_2DoF_AgentPositionMap(...
		x_all, AgentQuantity);
	Agent_ActualFormationObj.XData = plot_rx_set';
	Agent_ActualFormationObj.YData = plot_ry_set';

	[plot_rx_ref_set,plot_ry_ref_set] = Manipulator_2D_2DoF_AgentPositionMap(...
		kron(ones(AgentQuantity,1),xl_set(:,t_Nr)) + s_all_set(:,t_Nr), AgentQuantity);
	Agent_ReferenceFormationObj.XData = plot_rx_ref_set';
	Agent_ReferenceFormationObj.YData = plot_ry_ref_set';

	drawnow;
	pause(0.001);
end