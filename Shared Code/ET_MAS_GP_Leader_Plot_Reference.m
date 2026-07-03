%% Plot Referenve
% ET_MAS_GP_Leader_Plot_Reference
ReferenceFigureObj = figure('Name', 'Reference');
ReferenceStateAxesObj = subplot(2,2,[1,2],'Parent',ReferenceFigureObj);
plot(ReferenceStateAxesObj,x_l_set(1,:), x_l_set(2,:));
hold(ReferenceStateAxesObj, 'on');
ReferenceTrajectory_1_AxesObj = subplot(2,2,3,'Parent',ReferenceFigureObj);
plot(ReferenceTrajectory_1_AxesObj,t_set,x_l_set(1,:));
hold(ReferenceTrajectory_1_AxesObj,'on');
ReferenceTrajectory_2_AxesObj = subplot(2,2,4,'Parent',ReferenceFigureObj);
plot(ReferenceTrajectory_2_AxesObj,t_set,x_l_set(2,:));
hold(ReferenceTrajectory_2_AxesObj,'on');
for AgentNr = 1:AgentQuantity
	[s_i_set_plot, ~] = ET_MAS_GP_Leader_RelativePositionDynamics(t_set, AgentNr);
	
	plot(ReferenceStateAxesObj,x_l_set(1,:) + s_i_set_plot(1,:), x_l_set(2,:) + s_i_set_plot(2,:));
	plot(ReferenceTrajectory_1_AxesObj,t_set,x_l_set(1,:) + s_i_set_plot(1,:));
	plot(ReferenceTrajectory_2_AxesObj,t_set,x_l_set(2,:) + s_i_set_plot(2,:));
end
xlabel(ReferenceStateAxesObj,'x_1');
ylabel(ReferenceStateAxesObj,'x_2');
xlabel(ReferenceTrajectory_1_AxesObj,'t');
ylabel(ReferenceTrajectory_1_AxesObj,'x_1');
xlabel(ReferenceTrajectory_2_AxesObj,'t');
ylabel(ReferenceTrajectory_2_AxesObj,'x_2');
clear s_i_set_plot;