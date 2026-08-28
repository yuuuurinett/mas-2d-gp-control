function Summary = plot_ip_trigger_40s_summary(result_file)
%PLOT_IP_TRIGGER_40S_SUMMARY Per-second and cumulative packaged broadcasts.

data = load(result_file,'t_set','dac_inner_trigger_instance_set', ...
    'TrackingError_vector','prediction_error_norm_vector');
events = double(data.dac_inner_trigger_instance_set);
agent_count = size(events,1);
event_t = data.t_set(1:size(events,3));
duration = floor(data.t_set(end));
rate = nan(1,duration);
for second_i = 1:duration
    mask = event_t >= second_i-1 & event_t < second_i;
    rate(second_i) = sum(events(:,:,mask),'all')/agent_count;
end
centers = (0:duration-1)+0.5;
cumulative_per_agent = cumsum(sum(events,[1 2]))/agent_count;
cumulative_per_agent = reshape(cumulative_per_agent,1,[]);

MeanRate = mean(rate);
First5sRate = mean(rate(1:5));
Last5sRate = mean(rate(end-4:end));
Last10sRate = mean(rate(end-9:end));
TotalPerAgent = sum(events,'all')/agent_count;
MeanTrackingError = trapz(data.t_set,data.TrackingError_vector)/data.t_set(end);
MeanPredictionError = mean(data.prediction_error_norm_vector,'omitnan');
Summary = table(MeanRate,First5sRate,Last5sRate,Last10sRate, ...
    TotalPerAgent,MeanTrackingError,MeanPredictionError);

[folder,stem] = fileparts(result_file);
writetable(Summary,fullfile(folder,[stem,'_40s_summary.csv']));
fig = figure('Color','w','Position',[100 80 1120 700]);
layout = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax = nexttile(layout);
plot(ax,centers,rate,'LineWidth',1.6);
yline(ax,Last10sRate,'--','LineWidth',1.2, ...
    'DisplayName','Last 10 s mean');
xlabel(ax,'Time (s)'); ylabel(ax,'Broadcasts/agent/s');
title(ax,'Agent-level packaged broadcast rate'); grid(ax,'on'); box(ax,'on');
legend(ax,'Location','best');
ax = nexttile(layout);
plot(ax,event_t,cumulative_per_agent,'LineWidth',1.6);
xlabel(ax,'Time (s)'); ylabel(ax,'Cumulative broadcasts/agent');
title(ax,'Cumulative communication'); grid(ax,'on'); box(ax,'on');
title(layout,'M=600, R=1, epsilon=0.2, 40-s IP-DAC');
exportgraphics(fig,fullfile(folder,[stem,'_rate_summary.png']),'Resolution',230);
savefig(fig,fullfile(folder,[stem,'_rate_summary.fig']));
disp(Summary);
end
