function Summary = plot_ip_reachable_m600_comparison()
%PLOT_IP_REACHABLE_M600_COMPARISON M/epsilon accuracy-communication check.

repo_root = fileparts(fileparts(mfilename('fullpath')));
files = { ...
    fullfile(repo_root,'result','ip_reachable_uniform', ...
    'poe_M400_reachable_uniform_R1_eps_0p1_kappa_1_T30.mat'); ...
    fullfile(repo_root,'result','ip_reachable_m600', ...
    'poe_M600_reachable_R1_eps_0p1_kappa_1_T30.mat'); ...
    fullfile(repo_root,'result','ip_reachable_m600', ...
    'poe_M600_reachable_R1_eps_0p2_kappa_1_T30.mat')};
labels = {'M400, epsilon=0.1';'M600, epsilon=0.1';'M600, epsilon=0.2'};
M = [400;600;600]; Epsilon = [0.1;0.1;0.2];
n = numel(files);
MeanPredictionError = nan(n,1); MeanIPCENError = nan(n,1);
MeanTrackingError = nan(n,1); MeanBroadcastRate = nan(n,1);
First5sBroadcastRate = nan(n,1); Last5sBroadcastRate = nan(n,1);
data_set = cell(n,1); rates = cell(n,1);
for idx = 1:n
    data = load(files{idx}); data_set{idx} = data;
    valid = isfinite(data.prediction_error_norm_vector);
    MeanPredictionError(idx) = mean(data.prediction_error_norm_vector(valid));
    MeanIPCENError(idx) = mean(data.ip_cen_prediction_error_norm_vector(valid));
    MeanTrackingError(idx) = trapz( ...
        data.t_set,data.TrackingError_vector)/data.t_set(end);
    events = double(data.dac_inner_trigger_instance_set);
    event_t = data.t_set(1:size(events,3));
    rate = nan(1,30);
    for second_i = 1:30
        mask = event_t >= second_i-1 & event_t < second_i;
        rate(second_i) = sum(events(:,:,mask),'all')/size(events,1);
    end
    rates{idx} = rate;
    MeanBroadcastRate(idx) = mean(rate);
    First5sBroadcastRate(idx) = mean(rate(1:5));
    Last5sBroadcastRate(idx) = mean(rate(end-4:end));
end
Summary = table(labels,M,Epsilon,MeanPredictionError,MeanIPCENError, ...
    MeanTrackingError,MeanBroadcastRate,First5sBroadcastRate, ...
    Last5sBroadcastRate,'VariableNames',{'Case','M','Epsilon', ...
    'MeanIPDACError','MeanIPCENError','MeanTrackingError', ...
    'MeanBroadcastsPerAgentSecond','First5sBroadcastRate', ...
    'Last5sBroadcastRate'});
result_folder = fullfile(repo_root,'result','ip_reachable_m600');
writetable(Summary,fullfile(result_folder,'comparison.csv'));
save(fullfile(result_folder,'comparison.mat'),'Summary');

fig = figure('Color','w','Position',[80 70 1200 720]);
layout = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
colors = lines(n);
ax = nexttile(layout); hold(ax,'on');
for idx = 1:n
    plot(ax,data_set{idx}.t_set,data_set{idx}.TrackingError_vector, ...
        'Color',colors(idx,:),'LineWidth',1.5,'DisplayName',labels{idx});
end
xlabel(ax,'Time (s)'); ylabel(ax,'Tracking error norm');
title(ax,'Full-trajectory tracking'); grid(ax,'on'); box(ax,'on'); legend(ax);
ax = nexttile(layout); hold(ax,'on');
for idx = 1:n
    plot(ax,(0:29)+0.5,rates{idx},'Color',colors(idx,:), ...
        'LineWidth',1.5,'DisplayName',labels{idx});
end
xlabel(ax,'Time (s)'); ylabel(ax,'Broadcasts/agent/s');
title(ax,'Agent-level packaged broadcasts'); grid(ax,'on'); box(ax,'on');
legend(ax,'Location','best');
title(layout,'Reachable-domain IP-DAC: M and trigger threshold');
exportgraphics(fig,fullfile(result_folder, ...
    'm600_accuracy_trigger_comparison.png'),'Resolution',230);
savefig(fig,fullfile(result_folder,'m600_accuracy_trigger_comparison.fig'));
disp(Summary);
end
