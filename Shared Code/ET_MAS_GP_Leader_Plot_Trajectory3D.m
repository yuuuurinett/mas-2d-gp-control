%% Trajectory 3D
% ET_MAS_GP_Leader_Plot_Trajectory3D
FormationPlotInterval = 0.5;
%%
State3DFigureObj = figure('Name', 'Actual State 3D');
State3DAxesObj = subplot(1,2,1,'Parent',State3DFigureObj);
hold(State3DAxesObj,'on');
plot3(State3DAxesObj,x_l_set(1,:),x_l_set(2,:),t_set,'-.');
for AgentNr = 1:AgentQuantity
	x_i_set = x_set((AgentNr - 1) * x_dim + (1:x_dim),:);
	plot3(State3DAxesObj,x_i_set(1,:),x_i_set(2,:),t_set,'-');

	s_i_set = s_set((AgentNr - 1) * x_dim + (1:x_dim),:);
	r_i_set = x_l_set + s_i_set;
	plot3(State3DAxesObj,r_i_set(1,:),r_i_set(2,:),t_set,'--');
end
for FormationPlotTime = min(t_set):FormationPlotInterval:max(t_set)
	FormationPlotIndex = find(t_set >= FormationPlotTime,1);

	FormationPlot_x = x_set(:,FormationPlotIndex);
	FormationPlot_x = reshape(FormationPlot_x,x_dim,AgentQuantity);
	FormationPlot_x = [FormationPlot_x, FormationPlot_x(:,1)];
	plot3(State3DAxesObj,FormationPlot_x(1,:),FormationPlot_x(2,:), ...
		FormationPlotTime * ones(1,AgentQuantity + 1),'r.-');

	FormationPlot_x = s_set(:,FormationPlotIndex) + ...
		kron(ones(AgentQuantity,1), x_l_set(:,FormationPlotIndex));
	FormationPlot_x = reshape(FormationPlot_x,x_dim,AgentQuantity);
	FormationPlot_x = [FormationPlot_x, FormationPlot_x(:,1)];
	plot3(State3DAxesObj,FormationPlot_x(1,:),FormationPlot_x(2,:), ...
		FormationPlotTime * ones(1,AgentQuantity + 1),'b.--');
end
xlabel(State3DAxesObj,'x_1');
ylabel(State3DAxesObj,'x_2');
zlabel(State3DAxesObj,'t');
title(State3DAxesObj,'Actual Trajectory');
view(State3DAxesObj,[45,30]);
grid(State3DAxesObj,'on');
%%
Reference3DAxesObj = subplot(1,2,2,'Parent',State3DFigureObj);
hold(Reference3DAxesObj,'on');
plot3(Reference3DAxesObj,x_l_set(1,:),x_l_set(2,:),t_set,'-.');
for AgentNr = 1:AgentQuantity
	s_i_set = s_set((AgentNr - 1) * x_dim + (1:x_dim),:);
	r_i_set = x_l_set + s_i_set;
	plot3(Reference3DAxesObj,r_i_set(1,:),r_i_set(2,:),t_set,'--');
end
for FormationPlotTime = min(t_set):FormationPlotInterval:max(t_set)
	FormationPlotIndex = find(t_set >= FormationPlotTime,1);
	FormationPlot_x = s_set(:,FormationPlotIndex) + ...
		kron(ones(AgentQuantity,1), x_l_set(:,FormationPlotIndex));
	FormationPlot_x = reshape(FormationPlot_x,x_dim,AgentQuantity);
	FormationPlot_x = [FormationPlot_x, FormationPlot_x(:,1)];
	plot3(Reference3DAxesObj,FormationPlot_x(1,:),FormationPlot_x(2,:), ...
		FormationPlotTime * ones(1,AgentQuantity + 1),'.-');
end
xlabel(Reference3DAxesObj,'x_1');
ylabel(Reference3DAxesObj,'x_2');
zlabel(Reference3DAxesObj,'t');
title(Reference3DAxesObj,'Reference');
view(Reference3DAxesObj,[45,30]);
grid(Reference3DAxesObj,'on');
%%
clear r_i_set;
clear FormationPlotInterval FormationPlotTime FormationPlotIndex FormationPlot_x;