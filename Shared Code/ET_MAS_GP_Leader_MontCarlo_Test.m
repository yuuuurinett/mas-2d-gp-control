%% ET_MAS_GP_Leader_MontCarlo_Test
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
do_Randomization = true;
is_RandomInitialState = true;
MontCarloQuantity = 100;
if do_Simulation
	for EventTriggerTypeNr = [2,3]%1:numel(EventTriggerType_cell)
		mkdir(['Result\Event Trigger for MAS using GP with Leader\MontCarlo\', EventTriggerType_cell{EventTriggerTypeNr}])
		parfor MontCarloNr = 1:MontCarloQuantity
			SaveResultAddress = ['Result\Event Trigger for MAS using GP with Leader\MontCarlo\', ...
				EventTriggerType_cell{EventTriggerTypeNr},'\Test ',num2str(MontCarloNr),'.mat'];
			NotificationText = ['MontCarlo: <', EventTriggerType_cell{EventTriggerTypeNr},'-',num2str(MontCarloNr), '> finished!\n'];
			ET_MAS_GP_Leader_OnceComparison_TestFunc(EventTriggerType_cell{EventTriggerTypeNr}, ...
				do_Randomization,is_RandomInitialState, ...
				SaveResultAddress,NotificationText);
		end
	end
end
%% Result Analysis
close all;
%% Maximal Tracking Error
StartTime = 20;
MaximalTrackingErrorFigureObj = figure('Name','Maximal Tracking Error');
MaximalTrackingErrorAxesObj = axes(MaximalTrackingErrorFigureObj);
max_norm_vartheta_set = nan(numel(EventTriggerType_cell),MontCarloQuantity);
LegendSet = cell(numel(EventTriggerType_cell),1);
for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
	for MontCarloNr = 1:MontCarloQuantity
		SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\MontCarlo\', ...
				EventTriggerType_cell{EventTriggerTypeNr},'\Test ',num2str(MontCarloNr),'.mat']);
		vartheta_all_set = SimulationResult.vartheta_all_set;
		norm_vartheta_set = sqrt(sum(vartheta_all_set .^ 2));
		t_set = SimulationResult.t_set;

		norm_vartheta_set = norm_vartheta_set(t_set >= StartTime);
		max_norm_vartheta_set(EventTriggerTypeNr,MontCarloNr) = max(norm_vartheta_set);
	end
	LegendSet{EventTriggerTypeNr} = EventTriggerType_cell{EventTriggerTypeNr};
end
boxplot(MaximalTrackingErrorAxesObj,max_norm_vartheta_set');
set(MaximalTrackingErrorAxesObj,'YScale','log');
xticklabels(MaximalTrackingErrorAxesObj,LegendSet);
%
MethodNr_set = transpose(1:numel(EventTriggerType_cell));
tikz_variable_name_set = [];
for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
	eval([erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_max_e = ', ...
		'transpose(max_norm_vartheta_set(',num2str(EventTriggerTypeNr),',:));']);
	tikz_variable_name_set = [tikz_variable_name_set, ',', ...
		erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_max_e'];
end
data2txt_opt.fname = ['Result\Event Trigger for MAS using GP with Leader\MontCarlo\', ...
	'tikz format - txt file\','Method_max_e'];
eval(['data2txt(data2txt_opt,MethodNr_set',tikz_variable_name_set,');']);
% data2txt(data2txt_opt,MethodNr_set, ...
% 	NoTrigger_max_e,centralized_max_e,distributed_max_e,continuous_max_e,exact_max_e);
%
%% Number of Data Samples
NumberDataSetFigureObj = figure('Name','Number of Data Set');
NumberDataSetAxesObj = subplot(2,1,1,'Parent',NumberDataSetFigureObj);
TriggerAxesObj = subplot(2,1,2,'Parent',NumberDataSetFigureObj);
hold(NumberDataSetAxesObj,'on');
hold(TriggerAxesObj,'on');
LegendSet = cell(numel(EventTriggerType_cell),1);
for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
	for MontCarloNr = 1:MontCarloQuantity
		SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\MontCarlo\', ...
			EventTriggerType_cell{EventTriggerTypeNr},'\Test ',num2str(MontCarloNr),'.mat']);
		LocalGP_set = SimulationResult.LocalGP_set;
		Method_TriggerQuantity = sum(SimulationResult.trigger_set,2);
		for LocalGPNr = 1:numel(LocalGP_set)
			DataQuantitySet(LocalGPNr,EventTriggerTypeNr,MontCarloNr) = LocalGP_set{LocalGPNr}.DataQuantity;
			TriggerQuantitySet(:,EventTriggerTypeNr,MontCarloNr) = Method_TriggerQuantity;
		end
	end
	LegendSet{EventTriggerTypeNr} = EventTriggerType_cell{EventTriggerTypeNr};
end
mean_DataQuantitySet = mean(DataQuantitySet,3);
mean_TriggerQuantitySet = mean(TriggerQuantitySet,3);
variance_DataQuantitySet = std(DataQuantitySet,[],3);
variance_TriggerQuantitySet = std(TriggerQuantitySet,[],3);
DataQuantityBar = bar(NumberDataSetAxesObj,mean_DataQuantitySet);
TriggerQuantityBar = bar(TriggerAxesObj,mean_TriggerQuantitySet);
drawnow;
for DataQuantityBarNr = 1:numel(DataQuantityBar)
	ErrorBarXPosition = DataQuantityBar(DataQuantityBarNr).XData + DataQuantityBar(DataQuantityBarNr).XOffset;
	errorbar(NumberDataSetAxesObj,ErrorBarXPosition,mean_DataQuantitySet(:,DataQuantityBarNr),variance_DataQuantitySet(:,DataQuantityBarNr),'k.');
end
for TriggerQuantityBarNr = 1:numel(TriggerQuantityBar)
	ErrorBarXPosition = TriggerQuantityBar(TriggerQuantityBarNr).XData + TriggerQuantityBar(TriggerQuantityBarNr).XOffset;
	errorbar(TriggerAxesObj,ErrorBarXPosition,mean_TriggerQuantitySet(:,TriggerQuantityBarNr),variance_TriggerQuantitySet(:,TriggerQuantityBarNr),'k.');
end
title(NumberDataSetAxesObj, 'Size of Data Set');
title(TriggerAxesObj, 'Number of Triggers');
legend(NumberDataSetAxesObj,LegendSet,'NumColumns',numel(EventTriggerType_cell));
%
AgentNrSet = transpose(1:size(mean_DataQuantitySet,1));
tikz_variable_name_set_M = [];
tikz_variable_name_set_T = [];
for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
	Method_MaxDataSize_mean_set = mean_DataQuantitySet(:,EventTriggerTypeNr);
	Method_MaxDataSize_var_set  = variance_DataQuantitySet(:,EventTriggerTypeNr) ./ Method_MaxDataSize_mean_set;
	Method_MaxDataSize_min_set  = mean_DataQuantitySet(:,EventTriggerTypeNr) - variance_DataQuantitySet(:,EventTriggerTypeNr);
	Method_MaxDataSize_max_set  = mean_DataQuantitySet(:,EventTriggerTypeNr) + variance_DataQuantitySet(:,EventTriggerTypeNr);
	eval([erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_M_mean = Method_MaxDataSize_mean_set;']);
	eval([erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_M_var  = Method_MaxDataSize_var_set;']);
	eval([erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_M_min  = Method_MaxDataSize_min_set;']);
	eval([erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_M_max  = Method_MaxDataSize_max_set;']);
	tikz_variable_name_set_M = [tikz_variable_name_set_M,',', ...
		erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_M_mean,', ...
		erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_M_var,', ...
		erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_M_min,', ...
		erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_M_max'];

	Method_TriggerQuantity_mean_set = mean_TriggerQuantitySet(:,EventTriggerTypeNr);
	Method_TriggerQuantity_var_set  = variance_TriggerQuantitySet(:,EventTriggerTypeNr) ./ Method_TriggerQuantity_mean_set;
	Method_TriggerQuantity_min_set  = mean_TriggerQuantitySet(:,EventTriggerTypeNr) - variance_TriggerQuantitySet(:,EventTriggerTypeNr);
	Method_TriggerQuantity_max_set  = mean_TriggerQuantitySet(:,EventTriggerTypeNr) + variance_TriggerQuantitySet(:,EventTriggerTypeNr);
	eval([erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_T_mean = Method_TriggerQuantity_mean_set;']);
	eval([erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_T_var  = Method_TriggerQuantity_var_set;']);
	eval([erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_T_min  = Method_TriggerQuantity_min_set;']);
	eval([erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_T_max  = Method_TriggerQuantity_max_set;']);
	tikz_variable_name_set_T = [tikz_variable_name_set_T,',', ...
		erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_T_mean,', ...
		erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_T_var,', ...
		erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_T_min,', ...
		erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_T_max'];
end
data2txt_opt.fname = ['Result\Event Trigger for MAS using GP with Leader\MontCarlo\', ...
	'tikz format - txt file\','Method_M'];
eval(['data2txt(data2txt_opt,AgentNrSet',tikz_variable_name_set_M,tikz_variable_name_set_T,');']);
% data2txt(data2txt_opt,AgentNrSet, ...
% 	NoTrigger_M_mean,		NoTrigger_M_var,		NoTrigger_M_min,		NoTrigger_M_max, ...
% 	centralized_M_mean,		centralized_M_var,		centralized_M_min,		centralized_M_max, ...
% 	distributed_M_mean,		distributed_M_var,		distributed_M_min,		distributed_M_max, ...
% 	continuous001_M_mean,	continuous001_M_var,	continuous001_M_min,	continuous001_M_max, ...
% 	continuous0015_M_mean,	continuous0015_M_var,	continuous0015_M_min,	continuous0015_M_max, ...
% 	continuous002_M_mean,	continuous002_M_var,	continuous002_M_min,	continuous002_M_max, ...
% 	exact_M_mean,			exact_M_var,			exact_M_min,			exact_M_max, ...
% 	NoTrigger_T_mean,		NoTrigger_T_var,		NoTrigger_T_min,		NoTrigger_T_max, ...
% 	centralized_T_mean,		centralized_T_var,		centralized_T_min,		centralized_T_max, ...
% 	distributed_T_mean,		distributed_T_var,		distributed_T_min,		distributed_T_max, ...
% 	continuous001_T_mean,	continuous001_T_var,	continuous001_T_min,	continuous001_T_max, ...
% 	continuous0015_T_mean,	continuous0015_T_var,	continuous0015_T_min,	continuous0015_T_max, ...
% 	continuous002_T_mean,	continuous002_T_var,	continuous002_T_min,	continuous002_T_max, ...
% 	exact_T_mean,			exact_T_var,			exact_T_min,			exact_T_max);
%
% clear SimulationResult EventTriggerTypeNr LocalGP_set LocalGPNr DataQuantitySet Method_TriggerQuantity TriggerQuantitySet;
%% Overall Tracking Error over Time
StartTime = 20;
TrackingError_Time_FigureObj = figure('Name','Maximal Tracking Error');
TrackingError_Time_AxesObj = axes(TrackingError_Time_FigureObj);
LegendSet = cell(numel(EventTriggerType_cell),1);

SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\MontCarlo\', ...
	'distributed','\Test 1.mat']);
t_set_nominal = SimulationResult.t_set';
All_norm_vartheta_time_set = nan(MontCarloQuantity,numel(t_set_nominal),numel(EventTriggerType_cell));
Mean_norm_vartheta_time_set = nan(numel(EventTriggerType_cell),numel(t_set_nominal));
Std_norm_vartheta_time_set = nan(numel(EventTriggerType_cell),numel(t_set_nominal));
for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
	for MontCarloNr = 1:MontCarloQuantity
		SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\MontCarlo\', ...
				EventTriggerType_cell{EventTriggerTypeNr},'\Test ',num2str(MontCarloNr),'.mat']);
		vartheta_all_set = SimulationResult.vartheta_all_set;
		norm_vartheta_set = sqrt(sum(vartheta_all_set .^ 2));
		t_set = SimulationResult.t_set;
		norm_vartheta_nominal_set = interp1(t_set,norm_vartheta_set,t_set_nominal,'spline');
		All_norm_vartheta_time_set(MontCarloNr,:,EventTriggerTypeNr) = norm_vartheta_nominal_set;
	end
	Mean_norm_vartheta_time_set(EventTriggerTypeNr,:) = mean(All_norm_vartheta_time_set(:,:,EventTriggerTypeNr),1);
	Std_norm_vartheta_time_set(EventTriggerTypeNr,:) = std(All_norm_vartheta_time_set(:,:,EventTriggerTypeNr),1);
	LegendSet{EventTriggerTypeNr} = EventTriggerType_cell{EventTriggerTypeNr};
end
Max_norm_vartheta_time_set = Mean_norm_vartheta_time_set + Std_norm_vartheta_time_set;
Min_norm_vartheta_time_set = Mean_norm_vartheta_time_set - Std_norm_vartheta_time_set;
Min_norm_vartheta_time_set(Min_norm_vartheta_time_set < 1e-10) = 1e-10;

for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
	plot(TrackingError_Time_AxesObj,t_set_nominal,Mean_norm_vartheta_time_set(EventTriggerTypeNr,:), ...
		EventTriggerType_PlotConfiguration{EventTriggerTypeNr});
	hold(TrackingError_Time_AxesObj,'on');
end
for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
	fill(TrackingError_Time_AxesObj, ...
		[t_set_nominal;t_set_nominal(end:-1:1)], ...
		[Min_norm_vartheta_time_set(EventTriggerTypeNr,:),Max_norm_vartheta_time_set(EventTriggerTypeNr,end:-1:1)], ...
		EventTriggerType_PlotConfiguration{EventTriggerTypeNr}(1),'FaceAlpha',0.1,'EdgeAlpha',0);
end
legend(TrackingError_Time_AxesObj,LegendSet);
%
tikz_variable_name_set = [];
for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
	eval([erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_mean_e_t = ', ...
		'transpose(Mean_norm_vartheta_time_set(',num2str(EventTriggerTypeNr),',:));']);
	eval([erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_min_e_t = ', ...
		'transpose(Min_norm_vartheta_time_set(',num2str(EventTriggerTypeNr),',:));']);
	eval([erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_max_e_t = ', ...
		'transpose(Max_norm_vartheta_time_set(',num2str(EventTriggerTypeNr),',:));']);

	tikz_variable_name_set = [tikz_variable_name_set,',', ...
		erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_mean_e_t,', ...
		erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_min_e_t,', ...
		erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."]),'_max_e_t'];
end
data2txt_opt.fname = ['Result\Event Trigger for MAS using GP with Leader\MontCarlo\', ...
	'tikz format - txt file\','Method_e_over_t'];
eval(['data2txt(data2txt_opt',tikz_variable_name_set,',t_set_nominal);']);
%
%% Minimal Trigger Interval
% NumberDataSetFigureObj = figure('Name','Number of Data Set');
% NumberDataSetAxesObj = subplot(2,1,1,'Parent',NumberDataSetFigureObj);
% TriggerAxesObj = subplot(2,1,2,'Parent',NumberDataSetFigureObj);
% hold(NumberDataSetAxesObj,'on');
% hold(TriggerAxesObj,'on');

SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\MontCarlo\', ...
			EventTriggerType_cell{1},'\Test 1.mat']);
AgentQuantity = SimulationResult.AgentQuantity;

% LegendSet = cell(numel(EventTriggerType_cell),1);
min_t_trigger_interval_set = nan(AgentQuantity,MontCarloQuantity,numel(EventTriggerType_cell));
for EventTriggerTypeNr = [2,3]%1:numel(EventTriggerType_cell)
	for MontCarloNr = 1:MontCarloQuantity
		SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\MontCarlo\', ...
			EventTriggerType_cell{EventTriggerTypeNr},'\Test ',num2str(MontCarloNr),'.mat']);
		for AgentNr = 1:AgentQuantity
			t_trigger = find(SimulationResult.trigger_set(AgentNr,:) ~= 0) * SimulationResult.t_step;
			t_trigger_interval = t_trigger(2:end) - t_trigger(1:end-1);
			min_t_trigger_interval_set(AgentNr,MontCarloNr,EventTriggerTypeNr) = min([inf,t_trigger_interval]);
		end
	end
	LegendSet{EventTriggerTypeNr} = EventTriggerType_cell{EventTriggerTypeNr};
end
%% Single Agent Tracking Error
SingleAgentTrackingError_Time_FigureObj = figure('Name','Single Agent Tracking Error');
SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\MontCarlo\', ...
	'distributed','\Test 1.mat']);
AgentQuantity = SimulationResult.AgentQuantity;
t_set_nominal = SimulationResult.t_set';
x_dim = SimulationResult.x_dim;
SingleAgentTrackingError_Time_AxesSet = cell(AgentQuantity,1);
for AgentNr = 1:AgentQuantity
	SingleAgentTrackingError_Time_AxesSet{AgentNr} = subplot(AgentQuantity,1,AgentNr, ...
		'Parent',SingleAgentTrackingError_Time_FigureObj);
end
LegendSet = cell(numel(EventTriggerType_cell),1);

Agent_norm_vartheta_time_set = nan(MontCarloQuantity,numel(t_set_nominal),AgentQuantity,numel(EventTriggerType_cell));
Mean_Agent_norm_vartheta_time_set = nan(AgentQuantity,numel(t_set_nominal),numel(EventTriggerType_cell));
Std_Agent_norm_vartheta_time_set = nan(AgentQuantity,numel(t_set_nominal),numel(EventTriggerType_cell));
for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
	for MontCarloNr = 1:MontCarloQuantity
		SimulationResult = load(['Result\Event Trigger for MAS using GP with Leader\MontCarlo\', ...
				EventTriggerType_cell{EventTriggerTypeNr},'\Test ',num2str(MontCarloNr),'.mat']);
		t_set = SimulationResult.t_set;
		x_set = SimulationResult.x_set;
		s_set = SimulationResult.s_set;
		x_l_set = SimulationResult.x_l_set;
		e_set = s_set + repmat(x_l_set,AgentQuantity,1) - x_set;
		for AgentNr = 1:AgentQuantity
			e_i_set = e_set((AgentNr-1) * x_dim + (1:x_dim),:);
			norm_e_i_set = sqrt(sum(e_i_set .^ 2,1));
			norm_e_i_set = interp1(t_set,norm_e_i_set,t_set_nominal,'spline');
			Agent_norm_vartheta_time_set(MontCarloNr,:,AgentNr,EventTriggerTypeNr) = norm_e_i_set;
		end
	end
	for AgentNr = 1:AgentQuantity
		Mean_Agent_norm_vartheta_time_set(AgentNr,:,EventTriggerTypeNr) = ...
			mean(Agent_norm_vartheta_time_set(:,:,AgentNr,EventTriggerTypeNr),1);
		Std_Agent_norm_vartheta_time_set(AgentNr,:,EventTriggerTypeNr) = ...
			std(Agent_norm_vartheta_time_set(:,:,AgentNr,EventTriggerTypeNr),[],1);
	end
	LegendSet{EventTriggerTypeNr} = EventTriggerType_cell{EventTriggerTypeNr};
end
Max_Agent_norm_vartheta_time_set = Mean_Agent_norm_vartheta_time_set + Std_Agent_norm_vartheta_time_set;
Min_Agent_norm_vartheta_time_set = Mean_Agent_norm_vartheta_time_set - Std_Agent_norm_vartheta_time_set;
Min_Agent_norm_vartheta_time_set(Min_Agent_norm_vartheta_time_set < 1e-10) = 1e-10;
for AgentNr = 1:AgentQuantity
	for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
		semilogy(SingleAgentTrackingError_Time_AxesSet{AgentNr}, ...
			t_set_nominal,Mean_Agent_norm_vartheta_time_set(AgentNr,:,EventTriggerTypeNr), ...
			EventTriggerType_PlotConfiguration{EventTriggerTypeNr});
		hold(SingleAgentTrackingError_Time_AxesSet{AgentNr},'on');
	end
	for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
		fill(SingleAgentTrackingError_Time_AxesSet{AgentNr}, ...
			[t_set_nominal;t_set_nominal(end:-1:1)], ...
			[Min_Agent_norm_vartheta_time_set(AgentNr,:,EventTriggerTypeNr),Max_Agent_norm_vartheta_time_set(AgentNr,end:-1:1,EventTriggerTypeNr)], ...
			EventTriggerType_PlotConfiguration{EventTriggerTypeNr}(1),'FaceAlpha',0.1,'EdgeAlpha',0);
		hold(SingleAgentTrackingError_Time_AxesSet{AgentNr},'on');
	end
	legend(SingleAgentTrackingError_Time_AxesSet{AgentNr},LegendSet);
	ylim(SingleAgentTrackingError_Time_AxesSet{AgentNr},[1e-3,1e0]);
end
%

for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
	data2txt_opt.fname = ['Result\Event Trigger for MAS using GP with Leader\MontCarlo\', ...
		'tikz format - txt file\','AgentError_MonteCarlo_',erase(EventTriggerType_cell{EventTriggerTypeNr},[" ","."])];
	tikz_variable_name_set = [];
	for AgentNr = 1:AgentQuantity
		eval(['Mean_Agent',num2str(AgentNr), ...
			' = transpose(Mean_Agent_norm_vartheta_time_set(AgentNr,:,EventTriggerTypeNr));']);
		eval(['Min_Agent',num2str(AgentNr), ...
			' = transpose(Min_Agent_norm_vartheta_time_set(AgentNr,:,EventTriggerTypeNr));']);
		eval(['Max_Agent',num2str(AgentNr), ...
			' = transpose(Max_Agent_norm_vartheta_time_set(AgentNr,:,EventTriggerTypeNr));']);

		tikz_variable_name_set = [tikz_variable_name_set, ',', ...
			'Mean_Agent',num2str(AgentNr),',' ...
			'Min_Agent',num2str(AgentNr),',' ...
			'Max_Agent',num2str(AgentNr)];
	end
	eval(['data2txt(data2txt_opt',tikz_variable_name_set,',t_set_nominal);']);
end

% for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
% 	plot(TrackingError_Time_AxesObj,t_set_nominal,Mean_norm_vartheta_time_set(EventTriggerTypeNr,:), ...
% 		EventTriggerType_PlotConfiguration{EventTriggerTypeNr});
% 	hold(TrackingError_Time_AxesObj,'on');
% end
% for EventTriggerTypeNr = 1:numel(EventTriggerType_cell)
% 	fill(TrackingError_Time_AxesObj, ...
% 		[t_set_nominal;t_set_nominal(end:-1:1)], ...
% 		[Min_norm_vartheta_time_set(EventTriggerTypeNr,:),Max_norm_vartheta_time_set(EventTriggerTypeNr,end:-1:1)], ...
% 		EventTriggerType_PlotConfiguration{EventTriggerTypeNr}(1),'FaceAlpha',0.1,'EdgeAlpha',0);
% end
% legend(TrackingError_Time_AxesObj,LegendSet);






