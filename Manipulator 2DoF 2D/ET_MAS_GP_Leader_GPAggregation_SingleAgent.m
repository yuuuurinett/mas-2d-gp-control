function [f_hat_i, eta_aggregated_i,sigma_aggregated] = ET_MAS_GP_Leader_GPAggregation_SingleAgent( ...
	AgentNr, AgentBidirection_NeighbourSet, ...
	var_row_vector_SingleAgent, mu_row_cell_SingleAgent, beta, gamma)
%% GPOE
% var_vector_clean = rmmissing(var_row_vector_SingleAgent);
% var_i_aggregated = 1 / (sum(1 ./ var_vector_clean) / numel(var_vector_clean));
% AggregationWeight_set = var_row_vector_SingleAgent.^(-1) / var_i_aggregated^(-1) / numel(var_vector_clean);
% 
% 
% sigma_aggregated = rmmissing(sqrt(var_row_vector_SingleAgent)) * rmmissing(AggregationWeight_set');
% eta_aggregated_i = sqrt(beta) * sigma_aggregated + gamma;
% 
% f_hat_i = 0;
% for NeighborNr = 1:numel(var_row_vector_SingleAgent)
% 	if ~isnan(var_row_vector_SingleAgent(NeighborNr))
% 		f_hat_i = f_hat_i + AggregationWeight_set(NeighborNr) * mu_row_cell_SingleAgent{1, NeighborNr};
% 	end
% end
%% Smallest Error Bound
% [~,SelectedAgentIndex] = min(var_row_vector_SingleAgent);
% f_hat_i = mu_row_cell_SingleAgent{1, SelectedAgentIndex};
% sigma_aggregated = sqrt(var_row_vector_SingleAgent(SelectedAgentIndex));
% eta_aggregated_i = sqrt(beta) * sigma_aggregated + gamma;
%% Local
% f_hat_i = mu_row_cell_SingleAgent{1, AgentNr};
% sigma_aggregated = sqrt(var_row_vector_SingleAgent(AgentNr));
% eta_aggregated_i = sqrt(beta) * sigma_aggregated + gamma;
%% GPOE Posterior Variance
var_vector_clean = rmmissing(var_row_vector_SingleAgent);
var_i_aggregated = 1 / (sum(1 ./ var_vector_clean) / numel(var_vector_clean));
AggregationWeight_set = var_row_vector_SingleAgent.^(-1) / var_i_aggregated^(-1) / numel(var_vector_clean);


sigma_aggregated = sqrt(var_i_aggregated);
eta_aggregated_i = sqrt(beta) * sigma_aggregated + gamma;

f_hat_i = 0;
for NeighborNr = 1:numel(var_row_vector_SingleAgent)
	if ~isnan(var_row_vector_SingleAgent(NeighborNr))
		f_hat_i = f_hat_i + AggregationWeight_set(NeighborNr) * mu_row_cell_SingleAgent{1, NeighborNr};
	end
end
end