%% ET_MAS_GP_Leader_OnceComparison_Test
clc; clear;
%% Event Trigger Type
% NoTrigger, centralized, distributed, continuous
EventTriggerType_cell = {
	'NoTrigger';
	'centralized';
	'distributed';
	'continuous 0.01';
	'continuous 0.015';
	'continuous 0.02';
	'exact'
	};
EventTriggerType_PlotConfiguration = {
	'k--';
	'b-';
	'r-';
	'g-';
	'g--';
	'g.-';
	'm'
	};
%% Simulation
do_Simulation = true;
do_Randomization = false;
is_RandomInitialState = false;
if do_Simulation
	parfor EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
		SaveResultAddress = ['Result\Event Trigger for MAS using GP with Leader\Once Comparison\',EventTriggerType_cell{EventTriggerTypeNr},'.mat']
		NotificationText = ['Once Comparison: <', EventTriggerType_cell{EventTriggerTypeNr}, '> finished!\n']
		ET_MAS_GP_Leader_OnceComparison_TestFunc(EventTriggerType_cell{EventTriggerTypeNr}, ...
			do_Randomization,is_RandomInitialState, ...
			SaveResultAddress,NotificationText);
	end
end
%% Result Analysis
close all;
%% Tracking Error Curve
TrackingErrorFigureObj = figure('Name','Tracking Error');
TrackingErrorAxesObj = axes(TrackingErrorFigureObj);
for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
	SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
		EventTriggerType_cell{EventTriggerTypeNr},'.mat']);
	vartheta_all_set = SimulationResult.vartheta_all_set;
	norm_vartheta_set = sqrt(sum(vartheta_all_set .^ 2));
	t_set = SimulationResult.t_set;
	semilogy(TrackingErrorAxesObj,t_set,norm_vartheta_set);
	hold(TrackingErrorAxesObj,'on');

	eval(['norm_vartheta_set_tikz_',num2str(EventTriggerTypeNr),' = transpose(norm_vartheta_set);']);
end
legend(TrackingErrorAxesObj,EventTriggerType_cell,'NumColumns',4);
ylim(TrackingErrorAxesObj,[1e-3,1e0]);
xlabel(TrackingErrorAxesObj,'Time');
ylabel(TrackingErrorAxesObj,'Error');
%
t_set_tikz = t_set';
data2txt_opt.fname = ['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
	'tikz format - txt file\','Error'];
data2txt(data2txt_opt,t_set_tikz, ...
	norm_vartheta_set_tikz_1,norm_vartheta_set_tikz_2,norm_vartheta_set_tikz_3, ...
	norm_vartheta_set_tikz_4,norm_vartheta_set_tikz_5);
%
clear SimulationResult EventTriggerTypeNr vartheta_all_set norm_vartheta_set t_set;
%% Maximal Tracking Error
StartTime = 20;
MaximalTrackingErrorFigureObj = figure('Name','Maximal Tracking Error');
MaximalTrackingErrorAxesObj = axes(MaximalTrackingErrorFigureObj);
max_norm_vartheta_set = nan(numel(EventTriggerType_cell),1);
for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
	SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
		EventTriggerType_cell{EventTriggerTypeNr},'.mat']);
	vartheta_all_set = SimulationResult.vartheta_all_set;
	norm_vartheta_set = sqrt(sum(vartheta_all_set .^ 2));
	t_set = SimulationResult.t_set;

	norm_vartheta_set = norm_vartheta_set(t_set >= StartTime);
	max_norm_vartheta_set(EventTriggerTypeNr) = max(norm_vartheta_set);
end
bar(MaximalTrackingErrorAxesObj,max_norm_vartheta_set);
% set(MaximalTrackingErrorAxesObj,'YScale','log');
%% Number of Data Samples
NumberDataSetFigureObj = figure('Name','Number of Data Set');
NumberDataSetAxesObj = subplot(2,1,1,'Parent',NumberDataSetFigureObj);
TriggerAxesObj = subplot(2,1,2,'Parent',NumberDataSetFigureObj);
LegendSet = cell(numel(EventTriggerType_cell),1);
for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
	SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
		EventTriggerType_cell{EventTriggerTypeNr},'.mat']);
	LocalGP_set = SimulationResult.LocalGP_set;
	Method_TriggerQuantity = sum(SimulationResult.trigger_set,2);
	for LocalGPNr = 1:numel(LocalGP_set)
		DataQuantitySet(LocalGPNr,EventTriggerTypeNr) = LocalGP_set{LocalGPNr}.DataQuantity;
		TriggerQuantitySet(:,EventTriggerTypeNr) = Method_TriggerQuantity;
	end
	LegendSet{EventTriggerTypeNr} = EventTriggerType_cell{EventTriggerTypeNr};
end
bar(NumberDataSetAxesObj,DataQuantitySet);
bar(TriggerAxesObj,TriggerQuantitySet);
title(NumberDataSetAxesObj, 'Size of Data Set');
title(TriggerAxesObj, 'Number of Triggers');
legend(NumberDataSetAxesObj,LegendSet,'NumColumns',numel(EventTriggerType_cell));
clear SimulationResult EventTriggerTypeNr LocalGP_set LocalGPNr DataQuantitySet Method_TriggerQuantity TriggerQuantitySet;
%% Data Set
DataSetFigureObj = figure('Name','Data Set');
AxesRowQuantity = 2;
AxesColumnQuantity = 2;
AxesQuantity = AxesRowQuantity * AxesColumnQuantity;
DataSetAxesSet = cell(AxesRowQuantity * AxesColumnQuantity, 1);
for AxesRowNr = 1:AxesRowQuantity
	for AxesColumnNr = 1:AxesColumnQuantity
		DataSetAxesSet{(AxesRowNr - 1) * AxesColumnQuantity + AxesColumnNr} = ...
			subplot(AxesRowQuantity,AxesColumnQuantity, ...
			(AxesRowNr - 1) * AxesColumnQuantity + AxesColumnNr, ...
			'Parent',DataSetFigureObj);
		hold(DataSetAxesSet{(AxesRowNr - 1) * AxesColumnQuantity + AxesColumnNr}, 'on');
	end
end
MaximalLocalGPQuantity = AxesQuantity;
LegendSet = cell(numel(EventTriggerType_cell),1);
for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
	SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
		EventTriggerType_cell{EventTriggerTypeNr},'.mat']);
	LocalGP_set = SimulationResult.LocalGP_set;
	LocalGPQuantity = numel(LocalGP_set);
	for LocalGPNr = 1:numel(LocalGP_set)
		if LocalGPNr <= AxesQuantity
			plot(DataSetAxesSet{LocalGPNr}, ...
				LocalGP_set{LocalGPNr}.X(1,:), LocalGP_set{LocalGPNr}.X(2,:),'o');
		end
	end
	if LocalGPQuantity < MaximalLocalGPQuantity
		MaximalLocalGPQuantity = LocalGPQuantity;
	end
	LegendSet{EventTriggerTypeNr} = EventTriggerType_cell{EventTriggerTypeNr};
end
for MaximalLocalGPNr = 1:MaximalLocalGPQuantity
	title(DataSetAxesSet{MaximalLocalGPNr}, ...
			['Data Set ',num2str(MaximalLocalGPNr)]);
	legend(DataSetAxesSet{MaximalLocalGPNr}, LegendSet);
end
clear AxesRowQuantity AxesColumnQuantity AxesQuantity AxesRowNr AxesColumnNr SimulationResult;
clear LocalGP_set LocalGPNr LocalGPQuantity MaximalLocalGPQuantity MaximalLocalGPNr EventTriggerTypeNr;
%% Trigger Time
TriggerTimeFigureObj = figure('Name','Trigger Time');
TriggerTimeAxesObjSet = cell(numel(EventTriggerType_cell),1);
for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
	TriggerTimeAxesObjSet{EventTriggerTypeNr} = subplot(numel(EventTriggerType_cell),1,EventTriggerTypeNr, ...
		'Parent', TriggerTimeFigureObj);

	SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
		EventTriggerType_cell{EventTriggerTypeNr},'.mat']);
	LocalGP_set = SimulationResult.LocalGP_set;
	trigger_set = SimulationResult.trigger_set;
	t_set = SimulationResult.t_set;
	plot(TriggerTimeAxesObjSet{EventTriggerTypeNr}, ...
		t_set, 1 ./ trigger_set + [0:(numel(LocalGP_set) - 1)]', '*');
	xlim(TriggerTimeAxesObjSet{EventTriggerTypeNr},[min(t_set);max(t_set)]);
	%
	TriggerQuantitySet = sum(trigger_set,2);
	TriggerTimeSet = kron(ones(numel(LocalGP_set),1),t_set) .* (0 ./ trigger_set + 1);
	LegentSet = cell(numel(LocalGP_set),1);
	for LocalGPNr = 1:numel(LocalGP_set)
		TriggerTime = rmmissing(TriggerTimeSet(LocalGPNr,:));
% 		TriggerTime = TriggerTime(TriggerTime >= 0.01);
		TriggerInterval = TriggerTime(2:end) - TriggerTime(1:end-1);
		MinTriggerInterval = min(TriggerInterval);
		LegentSet{LocalGPNr} = ['Agent ', num2str(LocalGPNr), ...
			'(', num2str(TriggerQuantitySet(LocalGPNr)),'/',num2str(MinTriggerInterval), ')'];
		numel(find(TriggerInterval >= 0.015)) / numel(TriggerInterval) * 100
	end
	legend(TriggerTimeAxesObjSet{EventTriggerTypeNr}, LegentSet);
	title(TriggerTimeAxesObjSet{EventTriggerTypeNr}, EventTriggerType_cell{EventTriggerTypeNr});
	ylim(TriggerTimeAxesObjSet{EventTriggerTypeNr},[0.5; numel(LocalGP_set) + 0.5]);
end
%
for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
	SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
			EventTriggerType_cell{EventTriggerTypeNr},'.mat']);
	LocalGP_set = SimulationResult.LocalGP_set;
		trigger_set = SimulationResult.trigger_set;
		t_set = SimulationResult.t_set;
		for LocalGPNr = 1:numel(LocalGP_set)
			trigger_i = trigger_set(LocalGPNr,:)';
			time_i = t_set';

			time_i(trigger_i == 0) = [];
			trigger_i(trigger_i == 0) = [];
			trigger_i = trigger_i + LocalGPNr - 1;
			if strcmpi(EventTriggerType_cell{EventTriggerTypeNr}, 'centralized')
				eval(['trigger_centralized_',num2str(LocalGPNr),' = trigger_i;']);
				eval(['time_centralized_',num2str(LocalGPNr),' = time_i;']);
			elseif strcmpi(EventTriggerType_cell{EventTriggerTypeNr}, 'distributed')
				eval(['trigger_distributed_',num2str(LocalGPNr),' = trigger_i;']);
				eval(['time_distributed_',num2str(LocalGPNr),' = time_i;']);
			end
		end
end
%
data2txt_opt.fname = ['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
	'tikz format - txt file\','trigger'];
data2txt(data2txt_opt, ...
	time_centralized_1,trigger_centralized_1, ...
	time_centralized_2,trigger_centralized_2, ...
	time_centralized_3,trigger_centralized_3, ...
	time_centralized_4,trigger_centralized_4, ...
	time_distributed_1,trigger_distributed_1, ...
	time_distributed_2,trigger_distributed_2, ...
	time_distributed_3,trigger_distributed_3, ...
	time_distributed_4,trigger_distributed_4);
%% Reference
SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
		'NoTrigger','.mat']);
t_set = SimulationResult.t_set;
x_l_set = SimulationResult.x_l_set;
AgentQuantity = SimulationResult.AgentQuantity;
ET_MAS_GP_Leader_Plot_Reference;
% to txt
t_set_tikz = t_set';
x_l_1_set_tikz = x_l_set(1,:)';
x_l_2_set_tikz = x_l_set(2,:)';
for AgentNr = 1:AgentQuantity
	[s_i_set_plot, ~] = ET_MAS_GP_Leader_RelativePositionDynamics(t_set, AgentNr);
	eval(['s_',num2str(AgentNr),'_1_set_tikz = transpose(s_i_set_plot(1,:)) + x_l_1_set_tikz;']);
	eval(['s_',num2str(AgentNr),'_2_set_tikz = transpose(s_i_set_plot(2,:)) + x_l_2_set_tikz;']);
end
OfflineTrainingSamples = SimulationResult.LocalGP_set{1}.X;
dataset_x1 = OfflineTrainingSamples(1,:)';
dataset_x2 = OfflineTrainingSamples(2,:)';
data2txt_opt.fname = ['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
	'tikz format - txt file\','reference'];
data2txt_opt.ndata = 500;
data2txt(data2txt_opt,t_set_tikz,x_l_1_set_tikz,x_l_2_set_tikz, ...
	s_1_1_set_tikz,s_2_1_set_tikz,s_3_1_set_tikz,s_4_1_set_tikz, ...
	s_1_2_set_tikz,s_2_2_set_tikz,s_3_2_set_tikz,s_4_2_set_tikz, ...
	dataset_x1, dataset_x2);
% clear
clear SimulationResult t_set x_l_set AgentQuantity;
clear data2txt_opt;
%% 3D Trajectory
FormationPlotInterval = 2;
State3DFigureObj = figure('Name', 'Actual State 3D');
State3DAxesObjSet = cell(numel(EventTriggerType_cell),1);
for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
	State3DAxesObjSet{EventTriggerTypeNr} = subplot(...
		1,numel(EventTriggerType_cell),EventTriggerTypeNr, ...
		'Parent',State3DFigureObj);
	hold(State3DAxesObjSet{EventTriggerTypeNr},'on');

	SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
		EventTriggerType_cell{EventTriggerTypeNr},'.mat']);
	t_set = SimulationResult.t_set;
	s_set = SimulationResult.s_set;
	x_l_set = SimulationResult.x_l_set;
	x_set = SimulationResult.x_set;
	AgentQuantity = SimulationResult.AgentQuantity;
	x_dim = SimulationResult.x_dim;
	X_min = SimulationResult.X_min;
	X_max = SimulationResult.X_max;

	plot3(State3DAxesObjSet{EventTriggerTypeNr},x_l_set(1,:),x_l_set(2,:),t_set,'-.');
	for AgentNr = 1:AgentQuantity
		x_i_set = x_set((AgentNr - 1) * x_dim + (1:x_dim),:);
		plot3(State3DAxesObjSet{EventTriggerTypeNr},x_i_set(1,:),x_i_set(2,:),t_set,'-');

		s_i_set = s_set((AgentNr - 1) * x_dim + (1:x_dim),:);
		r_i_set = x_l_set + s_i_set;
		plot3(State3DAxesObjSet{EventTriggerTypeNr},r_i_set(1,:),r_i_set(2,:),t_set,'--');
	end
	for FormationPlotTime = min(t_set):FormationPlotInterval:max(t_set)
		FormationPlotIndex = find(t_set >= FormationPlotTime,1);

		FormationPlot_x = x_set(:,FormationPlotIndex);
		FormationPlot_x = reshape(FormationPlot_x,x_dim,AgentQuantity);
		FormationPlot_x = [FormationPlot_x, FormationPlot_x(:,1)];
		plot3(State3DAxesObjSet{EventTriggerTypeNr},FormationPlot_x(1,:),FormationPlot_x(2,:), ...
			FormationPlotTime * ones(1,AgentQuantity + 1),'r.-','LineWidth',1);

		FormationPlot_x = s_set(:,FormationPlotIndex) + ...
			kron(ones(AgentQuantity,1), x_l_set(:,FormationPlotIndex));
		FormationPlot_x = reshape(FormationPlot_x,x_dim,AgentQuantity);
		FormationPlot_x = [FormationPlot_x, FormationPlot_x(:,1)];
		plot3(State3DAxesObjSet{EventTriggerTypeNr},FormationPlot_x(1,:),FormationPlot_x(2,:), ...
			FormationPlotTime * ones(1,AgentQuantity + 1),'b.--','LineWidth',1);
	end
	xlabel(State3DAxesObjSet{EventTriggerTypeNr},'x_1');
	ylabel(State3DAxesObjSet{EventTriggerTypeNr},'x_2');
	zlabel(State3DAxesObjSet{EventTriggerTypeNr},'t');
	xlim(State3DAxesObjSet{EventTriggerTypeNr},[X_min(1);X_max(1)]);
	ylim(State3DAxesObjSet{EventTriggerTypeNr},[X_min(2);X_max(2)]);
	zlim(State3DAxesObjSet{EventTriggerTypeNr},[min(t_set);max(t_set)])
	title(State3DAxesObjSet{EventTriggerTypeNr},EventTriggerType_cell{EventTriggerTypeNr});
	view(State3DAxesObjSet{EventTriggerTypeNr},[45,30]);
	grid(State3DAxesObjSet{EventTriggerTypeNr},'on');
end


clear SimulationResult t_set s_set x_l_set x_set AgentQuantity x_dim;
clear AgentNr x_i_set s_i_set r_i_set;
clear FormationPlotTime FormationPlotIndex FormationPlot_x;
%% Topology

SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
		EventTriggerType_cell{1},'.mat']);
SimulationResult.MultiAgentSystem.Agent_Topology.show_topology;
%% Single Agent State and Error
AgentStateAndError_FigureSet = cell(numel(EventTriggerType_cell),1);
for EventTriggerTypeNr = [1,2,3,5,7]%1:numel(EventTriggerType_cell)
	EventTriggerTypeName = EventTriggerType_cell{EventTriggerTypeNr};
	AgentStateAndError_FigureSet{EventTriggerTypeNr}.FigureObj = figure('Name',EventTriggerTypeName);

	SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
		EventTriggerType_cell{EventTriggerTypeNr},'.mat']);
	t_set = SimulationResult.t_set';
	s_set = SimulationResult.s_set;
	x_l_set = SimulationResult.x_l_set';
	x_set = SimulationResult.x_set;
	AgentQuantity = SimulationResult.AgentQuantity;
	x_dim = SimulationResult.x_dim;
	X_min = SimulationResult.X_min;
	X_max = SimulationResult.X_max;

	AgentStateAndError_FigureSet{EventTriggerTypeNr}.AxesSet = cell(AgentQuantity,2);
	for AgentNr = 1:AgentQuantity
		x_i_set = x_set((AgentNr - 1) * x_dim + (1:x_dim),:)';
		s_i_set = s_set((AgentNr - 1) * x_dim + (1:x_dim),:)';
		r_i_set = x_l_set + s_i_set;
		e_i_set = x_i_set - r_i_set;
		norm_e_i_set = sqrt(sum(e_i_set .^ 2,2));

		AgentStateAndError_FigureSet{EventTriggerTypeNr}.AxesSet{AgentNr,1} = ...
			subplot(AgentQuantity,2,(AgentNr - 1) * 2 + 1, ...
			'Parent',AgentStateAndError_FigureSet{EventTriggerTypeNr}.FigureObj);
		plot(AgentStateAndError_FigureSet{EventTriggerTypeNr}.AxesSet{AgentNr,1}, ...
			t_set,x_i_set,'-',t_set,r_i_set,'--');

		AgentStateAndError_FigureSet{EventTriggerTypeNr}.AxesSet{AgentNr,2} = ...
			subplot(AgentQuantity,2,(AgentNr - 1) * 2 + 2, ...
			'Parent',AgentStateAndError_FigureSet{EventTriggerTypeNr}.FigureObj);
		semilogy(AgentStateAndError_FigureSet{EventTriggerTypeNr}.AxesSet{AgentNr,2}, ...
			t_set,abs(e_i_set),'-',t_set,norm_e_i_set,'--');
	end
end
%% Single Agent State and Error in One Figure
SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
		EventTriggerType_cell{1},'.mat']);
AgentQuantity = SimulationResult.AgentQuantity;
t_set_norminal = SimulationResult.t_set';
x_l_set = SimulationResult.x_l_set';
s_set = SimulationResult.s_set;
x_dim = SimulationResult.x_dim;
tikz_variable_name_set = [];
for AgentNr = 1:AgentQuantity
	s_i_set = s_set((AgentNr - 1) * x_dim + (1:x_dim),:);
	r_i_set = x_l_set' + s_i_set;
	for x_dim_Nr = 1:x_dim
		eval(['Agent',num2str(AgentNr),'_r',num2str(x_dim_Nr),' = transpose(r_i_set(x_dim_Nr,:));']);
		tikz_variable_name_set = [tikz_variable_name_set, ',', ...
			'Agent',num2str(AgentNr),'_r',num2str(x_dim_Nr)];
	end
end
x_l_1_set = x_l_set(:,1);
x_l_2_set = x_l_set(:,2);
data2txt_opt.fname = ['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
	'tikz format - txt file\','AgentStateError_Reference'];
eval(['data2txt(data2txt_opt,x_l_1_set,x_l_2_set',tikz_variable_name_set,',t_set_norminal);']);


AgentStateAndError_All_Figure = figure('Name','Agent State for All Methods');
AgentStateAndError_All_AxesSet = cell(AgentQuantity,2);
for AgentNr = 1:AgentQuantity
	s_i_set = s_set((AgentNr - 1) * x_dim + (1:x_dim),:)';
	r_i_set = x_l_set + s_i_set;

	AgentStateAndError_All_AxesSet{AgentNr,1} = ...
		subplot(AgentQuantity,2,(AgentNr - 1) * 2 + 1, ...
		'Parent',AgentStateAndError_All_Figure);
	plot(AgentStateAndError_All_AxesSet{AgentNr,1},t_set_norminal,r_i_set, ...
		'k--','LineWidth',1);
	hold(AgentStateAndError_All_AxesSet{AgentNr,1},'on');

	AgentStateAndError_All_AxesSet{AgentNr,2} = ...
		subplot(AgentQuantity,2,(AgentNr - 1) * 2 + 2, ...
		'Parent',AgentStateAndError_All_Figure);
end
for EventTriggerTypeNr = [1,2,3,5,7]%1:numel(EventTriggerType_cell)

	SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
		EventTriggerType_cell{EventTriggerTypeNr},'.mat']);
	t_set = SimulationResult.t_set;
	x_l_set = SimulationResult.x_l_set;
	s_set = SimulationResult.s_set;
	x_set = SimulationResult.x_set;
	
	x_l_set_norminal = nan(x_dim,numel(t_set_norminal));
% 	s_set_norminal = nan(x_dim * AgentQuantity,numel(t_set_norminal));
% 	x_set_norminal = nan(x_dim * AgentQuantity,numel(t_set_norminal));
	for x_dim_Nr = 1:x_dim
		x_l_set_norminal(x_dim_Nr,:) = interp1(x_l_set(x_dim_Nr,:),t_set,t_set_norminal,'spline');
% 		s_set_norminal = interp1(s_set,t_set,t_set_norminal,'spline');
% 		x_set_norminal = interp1(x_set,t_set,t_set_norminal,'spline');
	end

	data2txt_opt.fname = ['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
		'tikz format - txt file\',['AgentStateError_',erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."])]];
	tikz_variable_name_set = [];
	for AgentNr = 1:AgentQuantity
		x_i_set = x_set((AgentNr - 1) * x_dim + (1:x_dim),:);
		s_i_set = s_set((AgentNr - 1) * x_dim + (1:x_dim),:);
		r_i_set = x_l_set + s_i_set;
		e_i_set = x_i_set - r_i_set;
		norm_e_i_set = sqrt(sum(e_i_set .^ 2,1));

		x_i_set_norminal = nan(x_dim,numel(t_set_norminal));
		for x_dim_Nr = 1:x_dim
			x_i_set_norminal(x_dim_Nr,:) = interp1(t_set,x_i_set(x_dim_Nr,:),t_set_norminal','spline');
		end
		norm_e_i_set_norminal = interp1(t_set,norm_e_i_set,t_set_norminal,'spline');

		plot(AgentStateAndError_All_AxesSet{AgentNr,1}, ...
			t_set_norminal,x_i_set_norminal,[EventTriggerType_PlotConfiguration{EventTriggerTypeNr}(1),'-']);
		semilogy(AgentStateAndError_All_AxesSet{AgentNr,2}, ...
			t_set_norminal,norm_e_i_set_norminal,[EventTriggerType_PlotConfiguration{EventTriggerTypeNr}(1),'-']);
		hold(AgentStateAndError_All_AxesSet{AgentNr,2},'on');
	
		for x_dim_Nr = 1:x_dim
			eval(['Agent',num2str(AgentNr),'_x',num2str(x_dim_Nr),' = transpose(x_i_set_norminal(x_dim_Nr,:));']);
			tikz_variable_name_set = [tikz_variable_name_set, ',', ...
				'Agent',num2str(AgentNr),'_x',num2str(x_dim_Nr)];
		end
		eval(['Agent',num2str(AgentNr),'_e',' = transpose(norm_e_i_set);']);
		tikz_variable_name_set = [tikz_variable_name_set, ',', ...
				'Agent',num2str(AgentNr),'_e'];
		eval(['data2txt(data2txt_opt',tikz_variable_name_set,',t_set_norminal);']);

	end
end
%% Unknown Function
SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
		EventTriggerType_cell{1},'.mat']);
q_dim = SimulationResult.q_dim;
SystemOrder = SimulationResult.SystemOrder;
[X1_grid,X2_grid] = meshgrid(linspace(-1.5,1.5,80),linspace(-0.6,0.6,10));
X1_f = reshape(X1_grid,[],1);
X2_f = reshape(X2_grid,[],1);
Y_f = ET_MAS_GP_Leader_UnknownDynamics([X1_f,X2_f]',q_dim,SystemOrder)';
Y_grid = reshape(Y_f,10,80);
UnknownFunctionFigure = figure('Name','Unknown Function');
UnknownFunctionAxes = axes(UnknownFunctionFigure);
surf(UnknownFunctionAxes,X1_grid,X2_grid,Y_grid,'EdgeColor','none','FaceAlpha',0.5);
view(UnknownFunctionAxes,[-37.5,80]);
%
data2txt_opt.fname = ['Result\Event Trigger for MAS using GP with Leader\Once Comparison\', ...
		'tikz format - txt file\','UnknownFunction'];
data2txt_opt.ndata = 900;
data2txt(data2txt_opt, ...
	X1_f,X2_f,Y_f);

















