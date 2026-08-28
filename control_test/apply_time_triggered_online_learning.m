function [LocalGP_set, online_trigger_set, online_trigger_count] = ...
    apply_time_triggered_online_learning(LocalGP_set, online_trigger_set, ...
    online_trigger_count, t_Nr, x_all_matrix, unknown_scale, disturbance_scale)
% Add one local sample for every agent at a scheduled learning update.
% Each agent owns and counts its own GP dataset.  This helper implements the
% experiment's 0.1 s time-triggered acquisition, not the shared code's ET.

AgentQuantity = numel(LocalGP_set);
y_dim = LocalGP_set{1}.y_dim;

for AgentNr = 1:AgentQuantity
    x_i = x_all_matrix(:,AgentNr);
    y_i = unknown_scale * Manipulator_2D_2DoF_UnknownDynamics(x_i) + ...
        disturbance_scale * LocalGP_set{AgentNr}.SigmaN * randn(y_dim,1);

    if LocalGP_set{AgentNr}.DataQuantity >= LocalGP_set{AgentNr}.MaxDataQuantity
        LocalGP_set{AgentNr}.downdateParam(1);
    end

    LocalGP_set{AgentNr}.addPoint(x_i,y_i);
    online_trigger_set(AgentNr,t_Nr) = 1;
    online_trigger_count(AgentNr) = online_trigger_count(AgentNr) + 1;
end
end
