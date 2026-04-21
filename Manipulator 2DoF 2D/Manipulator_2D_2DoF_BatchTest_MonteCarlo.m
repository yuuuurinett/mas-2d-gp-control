%% Manipulator_2D_2DoF_BatchTest_MonteCarlo
clc; clear; close all;
%%
EventTriggerTypeSet = { ...
	'distributed'; 'centralized'; 'time'; 'offline'; 'exact'};
PlotColorSet = { ...
	'r'; 'b'; 'g'; 'c'; 'k'};
MonteCarloQuantity = 100;
SaveFolderParentName = ...
	['Result\Event Trigger for MAS using GP with Leader\Manipulator 2DoF 2D\', ...
	'Monte Carlo'];
%% Monte Carlo Test
do_simulation = false;
if do_simulation
	for EventTriggerTypeNr = 1:numel(EventTriggerTypeSet)
		EventTriggerType = EventTriggerTypeSet{EventTriggerTypeNr};
		SaveFolderName = [SaveFolderParentName,'\',EventTriggerType];
		parfor MonteCarloNr = 1:MonteCarloQuantity
			SaveFileName = num2str(MonteCarloNr);
			Manipulator_2D_2DoF_OnceComparison_TestFunc( ...
				EventTriggerType,SaveFolderName,SaveFileName);
			fprintf([EventTriggerType,' ', ...
				num2str(MonteCarloNr),'/',num2str(MonteCarloQuantity), ...
				' is finished! \n']);
		end
	end
end
%%
do_LoadResult = true;
if do_LoadResult
	StartTime = 20;
	RawResult = load([SaveFolderParentName,'\',EventTriggerTypeSet{1},'\1.mat'],'-mat', ...
		't_set','AgentQuantity');
	t_set = RawResult.t_set';
	AgentQuantity = RawResult.AgentQuantity;
	norm_vartheta_all_set_cell = cell(numel(EventTriggerTypeSet),1);
	max_norm_vartheta_all_set = nan(MonteCarloQuantity,numel(EventTriggerTypeSet));
	trigger_Quantity_set_cell = cell(numel(EventTriggerTypeSet),1);

	Mean_vartheta_all_set = nan(numel(t_set),numel(EventTriggerTypeSet));
	Std_vartheta_all_set = nan(numel(t_set),numel(EventTriggerTypeSet));
	Mean_trigger_Quantity_set = nan(AgentQuantity,numel(EventTriggerTypeSet));
	Std_trigger_Quantity_set = nan(AgentQuantity,numel(EventTriggerTypeSet));
	for EventTriggerTypeNr = 1:numel(EventTriggerTypeSet)
		EventTriggerType = EventTriggerTypeSet{EventTriggerTypeNr};
		SaveFolderName = [SaveFolderParentName,'\',EventTriggerType];
		norm_vartheta_all_set_cell{EventTriggerTypeNr} = ...
			nan(MonteCarloQuantity,numel(t_set));
		trigger_Quantity_set_cell{EventTriggerTypeNr} = ...
			nan(AgentQuantity,MonteCarloQuantity);
		for MonteCarloNr = 1:MonteCarloQuantity
			SaveFileName = num2str(MonteCarloNr);
			RawResult = load([SaveFolderName,'\',SaveFileName,'.mat'],'-mat', ...
				'vartheta_all_set','trigger_set','LocalGP_set');

			norm_vartheta_all_set = sqrt(sum(RawResult.vartheta_all_set .^ 2));
			norm_vartheta_all_set_cell{EventTriggerTypeNr}(MonteCarloNr,:) = ...
				norm_vartheta_all_set;
			max_norm_vartheta_all_set(MonteCarloNr,EventTriggerTypeNr) = ...
				max(norm_vartheta_all_set(t_set >= StartTime));
			
			if strcmpi(EventTriggerType,'offline')
				for AgentNr = 1:AgentQuantity
					trigger_Quantity_set_cell{EventTriggerTypeNr}(AgentNr,MonteCarloNr) = ...
						RawResult.LocalGP_set{AgentNr}.DataQuantity;
				end
			else
				trigger_Quantity_set_cell{EventTriggerTypeNr}(:,MonteCarloNr) = ...
					sum(RawResult.trigger_set,2);
			end

			fprintf([EventTriggerType,' ', ...
				num2str(MonteCarloNr),'/',num2str(MonteCarloQuantity), ...
				' is loaded! \n']);
		end
		Mean_vartheta_all_set(:,EventTriggerTypeNr) = mean(norm_vartheta_all_set_cell{EventTriggerTypeNr})';
		Std_vartheta_all_set(:,EventTriggerTypeNr) = std(norm_vartheta_all_set_cell{EventTriggerTypeNr})';

		Mean_trigger_Quantity_set(:,EventTriggerTypeNr) = mean(trigger_Quantity_set_cell{EventTriggerTypeNr},2);
		Std_trigger_Quantity_set(:,EventTriggerTypeNr) = std(trigger_Quantity_set_cell{EventTriggerTypeNr},[],2);
	end
	Min_vartheta_all_set = Mean_vartheta_all_set - Std_vartheta_all_set;
	Max_vartheta_all_set = Mean_vartheta_all_set + Std_vartheta_all_set;

	Min_trigger_Quantity_set = Mean_trigger_Quantity_set - Std_trigger_Quantity_set;
	Max_trigger_Quantity_set = Mean_trigger_Quantity_set + Std_trigger_Quantity_set;
	save([SaveFolderParentName,'\Result.mat'],'t_set','AgentQuantity', ...
		'norm_vartheta_all_set_cell','max_norm_vartheta_all_set','trigger_Quantity_set_cell', ...
		'Mean_vartheta_all_set','Std_vartheta_all_set','Min_vartheta_all_set','Max_vartheta_all_set', ...
		'Mean_trigger_Quantity_set','Std_trigger_Quantity_set','Min_trigger_Quantity_set','Max_trigger_Quantity_set');
else
	load([SaveFolderParentName,'\Result.mat']);
end
%% Tracking Error over Time
t_set = reshape(t_set,[],1);
FigureObj_OverallError = figure('Name','Overall Error');
AxesObj_OverallError = axes(FigureObj_OverallError);
for EventTriggerTypeNr = 1:numel(EventTriggerTypeSet)
	semilogy(AxesObj_OverallError,t_set,Mean_vartheta_all_set(:,EventTriggerTypeNr), ...
		'-','Color',PlotColorSet{EventTriggerTypeNr},'LineWidth',1);
	hold(AxesObj_OverallError,'on');
	fill(AxesObj_OverallError, ...
		[t_set;t_set(end:-1:1)], ...
		[Min_vartheta_all_set(:,EventTriggerTypeNr);Max_vartheta_all_set(end:-1:1,EventTriggerTypeNr)], ...
		PlotColorSet{EventTriggerTypeNr},'FaceAlpha',0.1,'EdgeAlpha',0);
end
set(AxesObj_OverallError,'YScale','log');
legend(AxesObj_OverallError,EventTriggerTypeSet,'NumColumns',5);
%
tikz_variable_name_set = [];
for EventTriggerTypeNr = 1:numel(EventTriggerTypeSet)
	eval([erase(EventTriggerTypeSet{EventTriggerTypeNr},[" ","."]),'_mean_e_t = ', ...
		'Mean_vartheta_all_set(:,',num2str(EventTriggerTypeNr),');']);
	eval([erase(EventTriggerTypeSet{EventTriggerTypeNr},[" ","."]),'_min_e_t = ', ...
		'Min_vartheta_all_set(:,',num2str(EventTriggerTypeNr),');']);
	eval([erase(EventTriggerTypeSet{EventTriggerTypeNr},[" ","."]),'_max_e_t = ', ...
		'Max_vartheta_all_set(:,',num2str(EventTriggerTypeNr),');']);

	tikz_variable_name_set = [tikz_variable_name_set,',', ...
		erase(EventTriggerTypeSet{EventTriggerTypeNr},[" ","."]),'_mean_e_t,', ...
		erase(EventTriggerTypeSet{EventTriggerTypeNr},[" ","."]),'_min_e_t,', ...
		erase(EventTriggerTypeSet{EventTriggerTypeNr},[" ","."]),'_max_e_t'];
end
data2txt_opt.fname = ['Result\Event Trigger for MAS using GP with Leader\Manipulator 2DoF 2D\Monte Carlo\', ...
	'tikz format - txt file\','Manipulator_Method_e_over_t'];
eval(['data2txt(data2txt_opt',tikz_variable_name_set,',t_set);']);

%% Max Error
FigureObj_MaximalTrackingError = figure('Name','Maximal Tracking Error');
AxesObj_MaximalTrackingError = axes(FigureObj_MaximalTrackingError);
boxplot(AxesObj_MaximalTrackingError,max_norm_vartheta_all_set);
set(AxesObj_MaximalTrackingError,'YScale','log');
xticklabels(AxesObj_MaximalTrackingError,EventTriggerTypeSet);
ylim(AxesObj_MaximalTrackingError,[0.07,1.3]);
%
MethodNr_set = transpose(1:numel(EventTriggerTypeSet));
tikz_variable_name_set = [];
for EventTriggerTypeNr = 1:numel(EventTriggerTypeSet)
	eval([erase(EventTriggerTypeSet{EventTriggerTypeNr},[" ","."]),'_max_e = ', ...
		'max_norm_vartheta_all_set(:,',num2str(EventTriggerTypeNr),');']);
	tikz_variable_name_set = [tikz_variable_name_set, ',', ...
		erase(EventTriggerTypeSet{EventTriggerTypeNr},[" ","."]),'_max_e'];
end
data2txt_opt.fname = ['Result\Event Trigger for MAS using GP with Leader\Manipulator 2DoF 2D\Monte Carlo\', ...
	'tikz format - txt file\','Manipulator_Method_max_e'];
eval(['data2txt(data2txt_opt,MethodNr_set',tikz_variable_name_set,');']);
%% Number of Trigger Events
FigureObj_NrTriggerEvent = figure('Name','Number of Trigger Events');
AxesObj_NrTriggerEvent = axes(FigureObj_NrTriggerEvent);
hold(AxesObj_NrTriggerEvent,'on');
BarObj_NrTriggerEvent = bar(AxesObj_NrTriggerEvent,Mean_trigger_Quantity_set);
ylim(AxesObj_NrTriggerEvent,[0,350]);
drawnow;
for EventTriggerTypeNr = 1:numel(EventTriggerTypeSet)
	ErrorBarXPosition = BarObj_NrTriggerEvent(EventTriggerTypeNr).XData + ...
		BarObj_NrTriggerEvent(EventTriggerTypeNr).XOffset;
	errorbar(AxesObj_NrTriggerEvent,ErrorBarXPosition, ...
		Mean_trigger_Quantity_set(:,EventTriggerTypeNr),Std_trigger_Quantity_set(:,EventTriggerTypeNr),'k.');
end
legend(AxesObj_NrTriggerEvent,EventTriggerTypeSet,'NumColumns',numel(EventTriggerTypeSet));
%
AgentNrSet = transpose(1:size(Mean_trigger_Quantity_set,1));
tikz_variable_name_set_M = [];
% tikz_variable_name_set_T = [];
for EventTriggerTypeNr = 1:numel(EventTriggerTypeSet)
	Method_MaxDataSize_mean_set = Mean_trigger_Quantity_set(:,EventTriggerTypeNr);
	Method_MaxDataSize_var_set  = Std_trigger_Quantity_set(:,EventTriggerTypeNr) ./ Method_MaxDataSize_mean_set;
	eval([erase(EventTriggerTypeSet{EventTriggerTypeNr},[" ","."]),'_M_mean = Method_MaxDataSize_mean_set;']);
	eval([erase(EventTriggerTypeSet{EventTriggerTypeNr},[" ","."]),'_M_var  = Method_MaxDataSize_var_set;']);

	tikz_variable_name_set_M = [tikz_variable_name_set_M,',', ...
		erase(EventTriggerTypeSet{EventTriggerTypeNr},[" ","."]),'_M_mean,', ...
		erase(EventTriggerTypeSet{EventTriggerTypeNr},[" ","."]),'_M_var'];
end
data2txt_opt.fname = ['Result\Event Trigger for MAS using GP with Leader\Manipulator 2DoF 2D\Monte Carlo\', ...
	'tikz format - txt file\','Manipulator_Method_M'];
eval(['data2txt(data2txt_opt,AgentNrSet',tikz_variable_name_set_M,');']);