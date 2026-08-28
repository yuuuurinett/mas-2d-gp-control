function Summary = plot_ip_agent_epsilon_sweep(epsilon_values)
%PLOT_IP_AGENT_EPSILON_SWEEP Compare packaged-trigger epsilon values.

if nargin < 1 || isempty(epsilon_values)
    epsilon_values = [0.005 0.02 0.04 0.08 0.1 0.2];
end
repo_root = fileparts(fileparts(mfilename('fullpath')));
result_folder = fullfile(repo_root,'result','ip_agent_epsilon_sweep');
if ~isfolder(result_folder), mkdir(result_folder); end

n = numel(epsilon_values);
data_set = cell(n,1);
for idx = 1:n
    epsilon_tag = strrep(sprintf('%.3g',epsilon_values(idx)),'.','p');
    if epsilon_values(idx) == 0.005
        source_file = fullfile(repo_root,'result','ip_dac_trigger_30s', ...
            'poe_random_M600_R10_T30.mat');
    else
        source_file = fullfile(result_folder,['poe_agent_eps_',epsilon_tag,'.mat']);
    end
    data_set{idx} = load(source_file);
end

Epsilon = epsilon_values(:);
MeanBroadcastsPerAgentSecond = nan(n,1);
EarlyBroadcastsPerAgentSecond = nan(n,1);
LateBroadcastsPerAgentSecond = nan(n,1);
MeanTrackingError = nan(n,1);
FinalTrackingError = nan(n,1);
MeanPredictionError = nan(n,1);
MeanDACTrackingError = nan(n,1);

colors = lines(n);
rate_fig = figure('Color','w','Position',[100 80 1120 680]);
layout = tiledlayout(rate_fig,2,1,'TileSpacing','compact','Padding','compact');
rate_ax = nexttile(layout); hold(rate_ax,'on');
tracking_ax = nexttile(layout); hold(tracking_ax,'on');

for idx = 1:n
    data = data_set{idx};
    events = double(data.dac_inner_trigger_instance_set);
    t = data.t_set(1:end-1);
    agent_count = size(events,1);
    window_count = ceil(data.t_set(end));
    rate = nan(1,window_count);
    for window_i = 1:window_count
        mask = t >= window_i-1 & t < window_i;
        rate(window_i) = sum(events(:,:,mask),'all')/agent_count;
    end
    centers = (0:window_count-1)+0.5;
    plot(rate_ax,centers,rate,'Color',colors(idx,:), ...
        'LineWidth',1.7,'DisplayName',sprintf('\\epsilon = %.3g',Epsilon(idx)));
    plot(tracking_ax,data.t_set,data.TrackingError_vector, ...
        'Color',colors(idx,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('\\epsilon = %.3g',Epsilon(idx)));

    duration = data.t_set(end)-data.t_set(1);
    MeanBroadcastsPerAgentSecond(idx) = sum(events,'all')/(agent_count*duration);
    EarlyBroadcastsPerAgentSecond(idx) = mean(rate(1:2));
    LateBroadcastsPerAgentSecond(idx) = mean(rate(end-1:end));
    MeanTrackingError(idx) = trapz(data.t_set,data.TrackingError_vector)/duration;
    FinalTrackingError(idx) = data.TrackingError_vector(end);
    MeanPredictionError(idx) = mean(data.prediction_error_norm_vector,'omitnan');
    MeanDACTrackingError(idx) = mean(data.dac_tracking_error_set,'omitnan');
end

xlabel(rate_ax,'Time (s)'); ylabel(rate_ax,'Broadcasts / agent / s');
title(rate_ax,'Agent-level packaged communication');
legend(rate_ax,'Location','best'); grid(rate_ax,'on'); box(rate_ax,'on');
xlabel(tracking_ax,'Time (s)'); ylabel(tracking_ax,'Tracking error');
title(tracking_ax,'Closed-loop tracking');
legend(tracking_ax,'Location','best'); grid(tracking_ax,'on'); box(tracking_ax,'on');
title(layout,'IP-DAC packaged-trigger threshold comparison');
exportgraphics(rate_fig,fullfile(result_folder,'epsilon_rate_tracking.png'),'Resolution',230);
savefig(rate_fig,fullfile(result_folder,'epsilon_rate_tracking.fig'));

%% Figure-3-style raster, one panel per epsilon.
raster_fig = figure('Color','w','Position',[100 50 1250 900]);
raster_layout = tiledlayout(raster_fig,n,1,'TileSpacing','compact','Padding','compact');
for idx = 1:n
    data = data_set{idx};
    events = logical(data.dac_inner_trigger_instance_set);
    t = data.t_set(1:size(events,3));
    ax = nexttile(raster_layout); hold(ax,'on');
    for agent_i = 1:size(events,1)
        [~,time_idx] = find(squeeze(events(agent_i,:,:)));
        plot(ax,t(time_idx),agent_i*ones(size(time_idx)), ...
            'k*','MarkerSize',2.2,'LineWidth',0.45,'HandleVisibility','off');
    end
    ylim(ax,[0.5,size(events,1)+0.5]); yticks(ax,1:size(events,1));
    ylabel(ax,'Agent'); xlim(ax,[0,t(end)]); grid(ax,'on'); box(ax,'on');
    title(ax,sprintf('\\epsilon = %.3g, mean %.1f broadcasts / agent / s', ...
        Epsilon(idx),MeanBroadcastsPerAgentSecond(idx)));
end
xlabel(raster_layout,'Time (s)');
title(raster_layout,'IP-DAC agent-level packaged broadcast times');
exportgraphics(raster_fig,fullfile(result_folder,'epsilon_broadcast_rasters.png'),'Resolution',230);
savefig(raster_fig,fullfile(result_folder,'epsilon_broadcast_rasters.fig'));

Summary = table(Epsilon,MeanBroadcastsPerAgentSecond, ...
    EarlyBroadcastsPerAgentSecond,LateBroadcastsPerAgentSecond, ...
    MeanTrackingError,FinalTrackingError,MeanPredictionError,MeanDACTrackingError);
writetable(Summary,fullfile(result_folder,'summary.csv'));
save(fullfile(result_folder,'summary.mat'),'Summary','epsilon_values');
disp(Summary);
end
