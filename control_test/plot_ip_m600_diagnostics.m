function plot_ip_m600_diagnostics(result_folder,m_value)
%PLOT_IP_M600_DIAGNOSTICS Trigger and dynamics diagnostics for IP AC/DAC.

if nargin < 1 || isempty(result_folder)
    repo_root = fileparts(fileparts(mfilename('fullpath')));
    result_folder = fullfile(repo_root,'result','ip_m_ablation');
end
if nargin < 2 || isempty(m_value), m_value = 600; end

dac = load(fullfile(result_folder,sprintf('ip_et_dac_poe_M%d.mat',m_value)));
ac = load(fullfile(result_folder,sprintf('ip_et_ac_poe_M%d.mat',m_value)));
agent_count = size(dac.f_hat_all_set,2);
colors = lines(agent_count);

%% Agent-level physical trigger-time raster
dac_events = squeeze(any(dac.dac_trigger_count_per_agent_point_set > 0,2));
if size(dac_events,1) ~= agent_count, dac_events = dac_events.'; end
ac_events = ac.ac_physical_trigger_count_set > 0;
if size(ac_events,1) ~= agent_count, ac_events = ac_events.'; end

trigger_fig = figure('Color','w','Position',[100 80 1200 720]);
layout = tiledlayout(trigger_fig,2,1,'TileSpacing','compact','Padding','compact');
plot_trigger_raster(nexttile(layout),dac.t_set(1:end-1),dac_events, ...
    colors,sprintf('IP-ET-DAC, M = %d',m_value));
plot_trigger_raster(nexttile(layout),ac.t_set(1:end-1),ac_events, ...
    colors,sprintf('IP-ET-AC, M = %d',m_value));
xlabel(layout,'Time (s)');
title(layout,'Agent-level communication activation times');
exportgraphics(trigger_fig,fullfile(result_folder, ...
    sprintf('trigger_time_raster_M%d.png',m_value)),'Resolution',220);
savefig(trigger_fig,fullfile(result_folder, ...
    sprintf('trigger_time_raster_M%d.fig',m_value)));

%% Continuous prediction-error comparison
prediction_fig = figure('Color','w','Position',[120 100 1100 620]);
ax = axes(prediction_fig); hold(ax,'on');
plot(ax,dac.t_set,dac.prediction_error_norm_vector,'LineWidth',1.8, ...
    'DisplayName','IP-ET-DAC');
plot(ax,ac.t_set,ac.prediction_error_norm_vector,'LineWidth',1.8, ...
    'DisplayName','IP-ET-AC');
xlabel(ax,'Time (s)');
ylabel(ax,'$\|f(x)-\hat{f}(x)\|_2$','Interpreter','latex');
title(ax,sprintf('Controller prediction error, M = %d',m_value));
legend(ax,'Location','best'); grid(ax,'on'); box(ax,'on');
xlim(ax,[dac.t_set(1),dac.t_set(end)]);
exportgraphics(prediction_fig,fullfile(result_folder, ...
    sprintf('prediction_error_M%d.png',m_value)),'Resolution',220);
savefig(prediction_fig,fullfile(result_folder, ...
    sprintf('prediction_error_M%d.fig',m_value)));

%% True dynamics versus predicted dynamics, six agents per panel
plot_dynamics_panels(dac,result_folder,m_value,'DAC',colors);
plot_dynamics_panels(ac,result_folder,m_value,'AC',colors);
end

function plot_trigger_raster(ax,t,event_set,colors,panel_title)
hold(ax,'on');
agent_count = size(event_set,1);
for agent_i = 1:agent_count
    event_times = t(logical(event_set(agent_i,:)));
    if ~isempty(event_times)
        scatter(ax,event_times,agent_i*ones(size(event_times)),18, ...
            colors(agent_i,:),'|','LineWidth',0.8);
    end
    text(ax,t(end)+0.08,agent_i,sprintf('n=%d',numel(event_times)), ...
        'Color',colors(agent_i,:),'FontSize',9,'VerticalAlignment','middle');
end
xlim(ax,[t(1),t(end)+0.75]); ylim(ax,[0.5,agent_count+0.5]);
yticks(ax,1:agent_count); ylabel(ax,'Agent'); title(ax,panel_title);
grid(ax,'on'); box(ax,'on');
end

function plot_dynamics_panels(data,result_folder,m_value,method_tag,colors)
fig = figure('Color','w','Position',[140 80 1250 760]);
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
y_dim = size(data.f_true_all_set,1);
agent_count = size(data.f_true_all_set,2);

for output_i = 1:y_dim
    true_ax = nexttile(layout,(output_i-1)*2+1); hold(true_ax,'on');
    pred_ax = nexttile(layout,(output_i-1)*2+2); hold(pred_ax,'on');
    for agent_i = 1:agent_count
        plot(true_ax,data.t_set,squeeze(data.f_true_all_set(output_i,agent_i,:)), ...
            'Color',colors(agent_i,:),'LineWidth',1.25, ...
            'DisplayName',sprintf('Agent %d',agent_i));
        plot(pred_ax,data.t_set,squeeze(data.f_hat_all_set(output_i,agent_i,:)), ...
            'Color',colors(agent_i,:),'LineWidth',1.25, ...
            'DisplayName',sprintf('Agent %d',agent_i));
    end
    title(true_ax,sprintf('True dynamics, output %d',output_i));
    title(pred_ax,sprintf('IP-%s prediction, output %d',method_tag,output_i));
    ylabel(true_ax,sprintf('f_%d',output_i));
    ylabel(pred_ax,sprintf('$\\hat{f}_{%d}$',output_i), ...
        'Interpreter','latex');
    grid(true_ax,'on'); box(true_ax,'on'); grid(pred_ax,'on'); box(pred_ax,'on');
    xlim(true_ax,[data.t_set(1),data.t_set(end)]);
    xlim(pred_ax,[data.t_set(1),data.t_set(end)]);
end
xlabel(layout,'Time (s)');
title(layout,sprintf('True and predicted unknown dynamics: IP-ET-%s, M = %d', ...
    method_tag,m_value));
legend(nexttile(layout,2),'Location','best');
exportgraphics(fig,fullfile(result_folder, ...
    sprintf('true_vs_prediction_%s_M%d.png',lower(method_tag),m_value)), ...
    'Resolution',220);
savefig(fig,fullfile(result_folder, ...
    sprintf('true_vs_prediction_%s_M%d.fig',lower(method_tag),m_value)));
end
