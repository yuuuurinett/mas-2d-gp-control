%% Manipulator_2D_2DoF_BatchTest_OneComparison
clc; clear; close all;
EventTriggerTypeSet = { ...
	'distributed'; 'centralized'; 'time'; 'offline'; 'exact'};
% EventTriggerTypeSet = { ...
% 	'distributed'; 'centralized'; 'offline'; 'exact'};
SaveFolderName = ...
	['Result\Event Trigger for MAS using GP with Leader\Manipulator 2DoF 2D\', ...
	'Once Comparison'];
%% Simulation
do_simulation = false;
if do_simulation
	for EventTriggerTypeNr = 1:numel(EventTriggerTypeSet)
		EventTriggerType = EventTriggerTypeSet{EventTriggerTypeNr};
		SaveFileName = EventTriggerType;
		rng(0);
		Manipulator_2D_2DoF_OnceComparison_TestFunc( ...
			EventTriggerType,SaveFolderName,SaveFileName);
	end
end
%% Load Saved Result
RawResult = load([SaveFolderName,'\',EventTriggerTypeSet{1},'.mat'],'-mat', ...
		'AgentQuantity');
AgentQuantity = RawResult.AgentQuantity;
%
norm_vartheta_set_cell = cell(numel(EventTriggerTypeSet),1);
norm_vartheta_all_set_cell = cell(numel(EventTriggerTypeSet),1);
x_set_cell = cell(numel(EventTriggerTypeSet),AgentQuantity);
trigger_set_cell = cell(numel(EventTriggerTypeSet),1);
for EventTriggerTypeNr = 1:numel(EventTriggerTypeSet)
	EventTriggerType = EventTriggerTypeSet{EventTriggerTypeNr};
	RawResult = load([SaveFolderName,'\',EventTriggerType,'.mat'],'-mat', ...
		'vartheta_all_set','trigger_set','x_dim','x_all_set');

% 	AgentQuantity = size(RawResult.vartheta_all_set,1) / RawResult.x_dim;
	norm_vartheta_set_cell{EventTriggerTypeNr} = nan(AgentQuantity,size(RawResult.vartheta_all_set,2));
	for AgentNr = 1:AgentQuantity
		norm_vartheta_set_cell{EventTriggerTypeNr}(AgentNr,:) = ...
			sqrt(sum(RawResult.vartheta_all_set((AgentNr - 1) * RawResult.x_dim + (1:RawResult.x_dim),:) .^ 2));
		x_set_cell{EventTriggerTypeNr,AgentNr} = ...
			RawResult.x_all_set((AgentNr - 1) * RawResult.x_dim + (1:RawResult.x_dim),:);
	end

	norm_vartheta_all_set_cell{EventTriggerTypeNr} = ...
		sqrt(sum(RawResult.vartheta_all_set .^ 2));
	trigger_set_cell{EventTriggerTypeNr} = RawResult.trigger_set;
end
%% Individual Error
FigureObj_IndividualError = figure('Name','Individual Error');
AxesObj_IndividualError_Set = cell(3,2);
RawResult = load([SaveFolderName,'\',EventTriggerTypeSet{1},'.mat'],'-mat', ...
		't_set','AgentQuantity');
t_set = RawResult.t_set';
AgentQuantity = RawResult.AgentQuantity;
for AgentNr = 1:AgentQuantity
	AxesObj_IndividualError_Set{AgentNr} = subplot(3,2,AgentNr,'Parent',FigureObj_IndividualError);
	for EventTriggerTypeNr = 1:numel(EventTriggerTypeSet)
		semilogy(AxesObj_IndividualError_Set{AgentNr},t_set,norm_vartheta_set_cell{EventTriggerTypeNr}(AgentNr,:),'-');
		hold(AxesObj_IndividualError_Set{AgentNr},'on');
	end
	title(AxesObj_IndividualError_Set{AgentNr},['Agent ',num2str(AgentNr)]);
end
%
for EventTriggerTypeNr = 1:numel(EventTriggerTypeSet)
	tikz_variable_name_set = [];
	data2txt_opt.fname = ['Result\Event Trigger for MAS using GP with Leader\Manipulator 2DoF 2D\Monte Carlo\', ...
		'tikz format - txt file\','Manipulator_AgentError_',EventTriggerTypeSet{EventTriggerTypeNr}];
	for AgentNr = 1:AgentQuantity
		eval(['Agent',num2str(AgentNr),'_e = transpose(norm_vartheta_set_cell{EventTriggerTypeNr}(AgentNr,:));']);
		tikz_variable_name_set = [tikz_variable_name_set, ',', ...
			'Agent',num2str(AgentNr),'_e'];
	end
	eval(['data2txt(data2txt_opt',tikz_variable_name_set,',t_set);']);
end
%% Overall Error
FigureObj_OverallError = figure('Name','Overall Error');
AxesObj_OverallError = axes(FigureObj_OverallError);

RawResult = load([SaveFolderName,'\',EventTriggerTypeSet{1},'.mat'],'-mat', ...
	't_set');
t_set = RawResult.t_set;
for EventTriggerTypeNr = 1:numel(EventTriggerTypeSet)
	semilogy(AxesObj_OverallError,t_set,norm_vartheta_all_set_cell{EventTriggerTypeNr},'-');
	hold(AxesObj_OverallError,'on');
end
legend(AxesObj_OverallError,EventTriggerTypeSet);
%% Trigger Events
FigureObj_TriggerEvent = figure('Name','Trigger Events');
% Centralized Event-trigger
AxesObj_CET_TriggerEvent = subplot(2,1,1,'Parent',FigureObj_TriggerEvent);
trigger_set_plot = trigger_set_cell{2};
AgentQuantity = size(trigger_set_plot,1);
trigger_set_plot(trigger_set_plot == 0) = nan;
sum(trigger_set_cell{2},2)
plot(AxesObj_CET_TriggerEvent, ...
	t_set,trigger_set_plot + kron(ones(1,numel(t_set)),(1:AgentQuantity)') - 1,'*');
axis(AxesObj_CET_TriggerEvent,[min(t_set),max(t_set),0.5,AgentQuantity + 0.5]);
legend(AxesObj_CET_TriggerEvent,cellstr(num2str(sum(trigger_set_cell{2},2))));
title(AxesObj_CET_TriggerEvent,'Centralized');
% Distributed Event-trigger
AxesObj_DET_TriggerEvent = subplot(2,1,2,'Parent',FigureObj_TriggerEvent);
trigger_set_plot = trigger_set_cell{1};
AgentQuantity = size(trigger_set_plot,1);
trigger_set_plot(trigger_set_plot == 0) = nan;
sum(trigger_set_cell{1},2)
plot(AxesObj_DET_TriggerEvent, ...
	t_set,trigger_set_plot + kron(ones(1,numel(t_set)),(1:AgentQuantity)') - 1,'*');
axis(AxesObj_DET_TriggerEvent,[min(t_set),max(t_set),0.5,AgentQuantity + 0.5]);
legend(AxesObj_DET_TriggerEvent,cellstr(num2str(sum(trigger_set_cell{1},2))));
title(AxesObj_DET_TriggerEvent,'Distributed');
%
tikz_variable_name_set = [];
RawResult = load([SaveFolderName,'\',EventTriggerTypeSet{1},'.mat'],'-mat', ...
		't_set','AgentQuantity');
t_set = RawResult.t_set';
AgentQuantity = RawResult.AgentQuantity;
for EventTriggerTypeNr = 1:numel(EventTriggerTypeSet)
	EventTriggerType = EventTriggerTypeSet{EventTriggerTypeNr};
	switch EventTriggerType
		case {'distributed','centralized'}
			trigger_set = trigger_set_cell{EventTriggerTypeNr};
			for AgentNr = 1:AgentQuantity
				time_i = t_set;
				trigger_i = trigger_set(AgentNr,:)';
				
				time_i(trigger_i == 0) = [];
				trigger_i(trigger_i == 0) = [];
				trigger_i = trigger_i + AgentNr - 1;

				eval(['trigger_',EventTriggerType,'_',num2str(AgentNr),' = trigger_i;']);
				eval(['time_',EventTriggerType,'_',num2str(AgentNr),' = time_i;']);

				tikz_variable_name_set = [tikz_variable_name_set, ',', ...
					'trigger_',EventTriggerType,'_',num2str(AgentNr),',', ...
					'time_',EventTriggerType,'_',num2str(AgentNr)];
			end
	end
end
data2txt_opt.fname = ['Result\Event Trigger for MAS using GP with Leader\Manipulator 2DoF 2D\Monte Carlo\', ...
	'tikz format - txt file\','Manipulator_trigger'];
eval(['data2txt(data2txt_opt',tikz_variable_name_set,');']);
% data2txt(data2txt_opt, ...
% 	time_centralized_1,trigger_centralized_1, ...
% 	time_centralized_2,trigger_centralized_2, ...
% 	time_centralized_3,trigger_centralized_3, ...
% 	time_centralized_4,trigger_centralized_4, ...
% 	time_distributed_1,trigger_distributed_1, ...
% 	time_distributed_2,trigger_distributed_2, ...
% 	time_distributed_3,trigger_distributed_3, ...
% 	time_distributed_4,trigger_distributed_4);