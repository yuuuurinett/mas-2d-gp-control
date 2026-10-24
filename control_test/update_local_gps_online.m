function [LocalGP_set, update_history, update_count, eviction_count] = ...
    update_local_gps_online(LocalGP_set, update_history, update_count, ...
    sample_index, x_all_matrix, unknown_scale, disturbance_scale, update_mask)
%UPDATE_LOCAL_GPS_ONLINE Add one owned sample to selected local GP datasets.
%
% This is the only online data-acquisition path shared by LOCAL, CEN, NBR,
% TP and IP.  Scheduling and communication are intentionally kept outside:
% callers decide *when* an update is due and, for event-triggered methods,
% which agents are selected through update_mask.
%
% FIFO eviction (below) only fires once an agent's LocalGP already holds
% MaxDataQuantity points and a NEW point is about to be added, i.e. once
% that agent's update count so far exceeds MaxDataQuantity. With the
% current MC15 settings (online update every 0.1 s over a 20 s horizon =>
% 200 updates/agent, MaxDataQuantity == 200) that condition is never met,
% so eviction_count is always 0 and this branch is a no-op for this batch.
%
% Do not delete it: it exists for configurations where the update budget
% DOES exceed MaxDataQuantity (e.g. an M-ablation run with a smaller cap,
% a longer horizon, or a shorter update interval). Without it, addPoint
% would simply refuse the new point once full (flag=-1) instead of
% keeping a rolling window of the most recent data.
%
% If eviction_count is ever nonzero on a run where you expected it to
% stay within budget, that is a signal the cap/horizon/update-interval
% assumption changed upstream and should be checked before trusting the
% results -- see the warning below.

agent_count = numel(LocalGP_set);
if nargin < 8 || isempty(update_mask)
    update_mask = true(agent_count,1);
else
    update_mask = logical(update_mask(:));
end
if numel(update_mask) ~= agent_count
    error('update_mask must contain one logical value per agent.');
end

eviction_count = zeros(agent_count,1);
for agent_idx = find(update_mask).'
    gp = LocalGP_set{agent_idx};

    new_x = x_all_matrix(:,agent_idx);
    new_y = unknown_scale * Manipulator_2D_2DoF_UnknownDynamics(new_x) + ...
        disturbance_scale * gp.SigmaN * randn(gp.y_dim,1);

    is_full = gp.DataQuantity >= gp.MaxDataQuantity;
    if is_full
        % FIFO: drop the oldest point so this new one fits within the
        % fixed per-agent data budget (MaxDataQuantity). Under the
        % current MC15 settings this should never execute -- see note
        % above.
        gp.downdateParam(1);
        eviction_count(agent_idx) = 1;
    end
    gp.addPoint(new_x,new_y);

    update_history(agent_idx,sample_index) = 1;
    update_count(agent_idx) = update_count(agent_idx) + 1;
end

if any(eviction_count)
    warning('update_local_gps_online:unexpectedEviction', ...
        ['FIFO eviction triggered for %d agent(s) at sample %d. If ' ...
        'every agent''s update budget was expected to stay within ' ...
        'MaxDataQuantity for this run, check MaxDataQuantity, the ' ...
        'simulation horizon, and the online-update interval -- this is ' ...
        'the same class of silent-mismatch issue that produced the ' ...
        'earlier MC15_BIG_TABLE_formation discrepancy.'], ...
        nnz(eviction_count),sample_index);
end
end
