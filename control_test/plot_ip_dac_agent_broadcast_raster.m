function Counts = plot_ip_dac_agent_broadcast_raster(result_file)
%PLOT_IP_DAC_AGENT_BROADCAST_RASTER Figure-3-style communication raster.
%
% One asterisk denotes one agent-level packaged broadcast in one DAC inner
% round. Multiple inner-round broadcasts at the same physical time overlap;
% the exact per-agent total is printed to the right of the axes.

if nargin < 1 || isempty(result_file)
    repo_root = fileparts(fileparts(mfilename('fullpath')));
    result_file = fullfile(repo_root,'result','ip_dac_trigger_30s', ...
        'poe_random_M600_R10_T30.mat');
end
data = load(result_file,'t_set','dac_inner_trigger_instance_set', ...
    'DACStepsPerPhysicalStep','NumInducingPoints','DACTriggerEpsilon');
[result_folder,result_stem] = fileparts(result_file);

events = logical(data.dac_inner_trigger_instance_set);
[agent_count,round_count,step_count] = size(events);
t = data.t_set(1:step_count);
Counts = squeeze(sum(events,[2 3]));

fig = figure('Color','w','Position',[100 100 1220 500]);
ax = axes(fig); hold(ax,'on');
agent_colors = lines(agent_count);
for agent_i = 1:agent_count
    % Keep the [DAC round x physical step] orientation even when there is
    % only one DAC round; squeeze() would turn a 1-by-K slice into K-by-1.
    agent_events = reshape(events(agent_i,:,:),round_count,step_count);
    [~,time_idx] = find(agent_events);
    event_times = t(time_idx);
    plot(ax,event_times,agent_i*ones(size(event_times)), ...
        '*','Color',agent_colors(agent_i,:), ...
        'MarkerSize',3.2,'LineWidth',0.55, ...
        'HandleVisibility','off');
    text(ax,t(end)+0.35,agent_i,sprintf('n = %d',Counts(agent_i)), ...
        'Color',agent_colors(agent_i,:),'FontSize',10, ...
        'VerticalAlignment','middle');
end

xlim(ax,[t(1),t(end)+1.7]);
ylim(ax,[0.5,agent_count+0.5]);
yticks(ax,1:agent_count);
yticklabels(ax,compose('Agent %d',1:agent_count));
xlabel(ax,'Time (s)');
ylabel(ax,'Agent');
title(ax,sprintf(['IP-DAC agent-level packaged broadcast times ', ...
    '(M=%d, rounds=%d, \\epsilon=%.3g)'], ...
    data.NumInducingPoints,round_count,data.DACTriggerEpsilon));
grid(ax,'on'); box(ax,'on');
set(ax,'GridAlpha',0.14,'MinorGridAlpha',0.08);

output_base = fullfile(result_folder,[result_stem,'_agent_broadcast_raster']);
exportgraphics(fig,[output_base,'.png'],'Resolution',240);
savefig(fig,[output_base,'.fig']);
writetable(table((1:agent_count)',Counts, ...
    'VariableNames',{'Agent','PackagedBroadcastCount'}), ...
    [output_base,'.csv']);
end
