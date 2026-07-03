%% State
% ET_MAS_GP_Leader_Plot_SimulationState
StateFigureObj = figure('Name', 'Actual State');
%
StatePhaseAxesObj = subplot(2,2,[1,2],'Parent',StateFigureObj);
hold(StatePhaseAxesObj,'on');
plot(StatePhaseAxesObj,x_l_set(1,:), x_l_set(2,:),'-.');
Trajectory_1_AxesObj = subplot(2,2,3,'Parent',StateFigureObj);
hold(Trajectory_1_AxesObj, 'on');
plot(Trajectory_1_AxesObj,t_set, x_l_set(1,:),'-.');
Trajectory_2_AxesObj = subplot(2,2,4,'Parent',StateFigureObj);
hold(Trajectory_2_AxesObj, 'on');
plot(Trajectory_2_AxesObj,t_set, x_l_set(2,:),'-.');
StateLegendSet = {['Leader']};
for AgentNr = 1:AgentQuantity
	x_i_set = x_set((AgentNr - 1) * x_dim + (1:x_dim),:);
	plot(StatePhaseAxesObj,x_i_set(1,:), x_i_set(2,:),'-');
	plot(Trajectory_1_AxesObj,t_set, x_i_set(1,:));
	plot(Trajectory_2_AxesObj,t_set, x_i_set(2,:));
	StateLegendSet{numel(StateLegendSet) + 1} = ['x_',num2str(AgentNr)];

	x_r_i_set = s_set((AgentNr - 1) * x_dim + (1:x_dim),:) + x_l_set;
	plot(StatePhaseAxesObj,x_r_i_set(1,:), x_r_i_set(2,:),'--');
	plot(Trajectory_1_AxesObj,t_set,x_r_i_set(1,:),'--');
	plot(Trajectory_2_AxesObj,t_set,x_r_i_set(2,:),'--');
	StateLegendSet{numel(StateLegendSet) + 1} = ['s_',num2str(AgentNr),' + x_l'];
end
legend(StatePhaseAxesObj,StateLegendSet);
legend(Trajectory_1_AxesObj,StateLegendSet);
legend(Trajectory_2_AxesObj,StateLegendSet);