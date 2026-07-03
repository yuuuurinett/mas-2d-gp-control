%% Trigger
% ET_MAS_GP_Leader_Plot_Trigger
TriggerFigureObj = figure('Name', 'Trigger');
TriggerAxesObj = axes(TriggerFigureObj);
plot(TriggerAxesObj,t_set,trigger_set + (0:(AgentQuantity-1))');
TriggerQuantity = sum(trigger_set,2);
TriggerLegentSet = cell(AgentQuantity,1);
for AgentNr = 1:AgentQuantity
	TriggerLegentSet{AgentNr} = ['Trigger for Agent ', num2str(AgentNr), ' (',num2str(TriggerQuantity(AgentNr)),')'];
end
legend(TriggerAxesObj, TriggerLegentSet);
title(TriggerAxesObj,EventTriggerType);