function plot_control_task_requirements(ResultFolder)
% Requirement-aligned control plots. Multi-output dynamics are shown by
% vector norms rather than separate output-dimension panels.

if nargin < 1 || isempty(ResultFolder)
    ProjectRoot = fileparts(fileparts(mfilename('fullpath')));
    ResultFolder = fullfile(ProjectRoot,'result','tp_ip_every_step_baseline');
end

TPF = load(fullfile(ResultFolder,'test_point_poe_every_step.mat'));
EXF = load(fullfile(ResultFolder,'exact_true_dynamics.mat'));
TPN = load(fullfile(ResultFolder,'test_point_poe_every_step_noformation.mat'));
EXN = load(fullfile(ResultFolder,'exact_true_dynamics_noformation.mat'));

%% Tracking error: with and without formation
fig = figure('Color','w','Position',[60,60,1250,540]);
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
plot_tracking_panel(nexttile(layout),TPF,EXF,'With formation');
plot_tracking_panel(nexttile(layout),TPN,EXN,'Without formation');
title(layout,'Test-point online GP tracking performance');
export_clean(fig);
exportgraphics(fig,fullfile(ResultFolder, ...
    'control_tracking_with_without_formation.png'),'Resolution',220);
savefig(fig,fullfile(ResultFolder, ...
    'control_tracking_with_without_formation.fig'));

%% TP prediction versus true dynamics: norm across the two output dimensions
fig = figure('Color','w','Position',[80,50,1250,760]);
layout = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(layout); hold(ax1,'on');
ax2 = nexttile(layout); hold(ax2,'on');
colors = lines(size(TPF.f_true_all_set,2));
for agent_i = 1:size(TPF.f_true_all_set,2)
    true_norm = squeeze(vecnorm(TPF.f_true_all_set(:,agent_i,:),2,1));
    pred_norm = squeeze(vecnorm(TPF.f_hat_all_set(:,agent_i,:),2,1));
    error_norm = squeeze(vecnorm( ...
        TPF.f_hat_all_set(:,agent_i,:)-TPF.f_true_all_set(:,agent_i,:),2,1));
    plot(ax1,TPF.t_set,true_norm,'-','Color',colors(agent_i,:), ...
        'LineWidth',1.1,'DisplayName',sprintf('Agent %d true',agent_i));
    plot(ax1,TPF.t_set,pred_norm,'--','Color',colors(agent_i,:), ...
        'LineWidth',1.1,'DisplayName',sprintf('Agent %d prediction',agent_i));
    semilogy(ax2,TPF.t_set,max(error_norm,eps),'-', ...
        'Color',colors(agent_i,:),'LineWidth',1.1, ...
        'DisplayName',sprintf('Agent %d',agent_i));
end
grid(ax1,'on'); box(ax1,'on');
xlabel(ax1,'Time (s)'); ylabel(ax1,'Dynamics vector norm');
title(ax1,'True and predicted unknown-dynamics magnitude');
legend(ax1,'Location','eastoutside');
grid(ax2,'on'); box(ax2,'on');
xlabel(ax2,'Time (s)'); ylabel(ax2,'Prediction-error norm');
title(ax2,'Per-agent vector prediction error');
legend(ax2,'Location','eastoutside');
title(layout,'Test-point prediction versus true dynamics (two outputs combined)');
export_clean(fig);
exportgraphics(fig,fullfile(ResultFolder, ...
    'control_tp_prediction_vs_true_vector_norm.png'),'Resolution',220);
savefig(fig,fullfile(ResultFolder, ...
    'control_tp_prediction_vs_true_vector_norm.fig'));

%% Prediction error and tracking error on a shared time axis
fig = figure('Color','w','Position',[80,80,1100,650]);
layout = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(layout);
semilogy(ax1,TPF.t_set,max(TPF.prediction_error_norm_vector,eps), ...
    'LineWidth',1.6);
grid(ax1,'on'); box(ax1,'on'); ylabel(ax1,'Prediction-error norm');
title(ax1,'Aggregated unknown-dynamics prediction error');
ax2 = nexttile(layout);
semilogy(ax2,TPF.t_set,max(TPF.TrackingError_vector,eps), ...
    'LineWidth',1.6);
grid(ax2,'on'); box(ax2,'on');
xlabel(ax2,'Time (s)'); ylabel(ax2,'Tracking error ||\vartheta||');
title(ax2,'Closed-loop tracking error');
title(layout,'Prediction error versus tracking error');
export_clean(fig);
exportgraphics(fig,fullfile(ResultFolder, ...
    'control_prediction_error_vs_tracking_error.png'),'Resolution',220);
savefig(fig,fullfile(ResultFolder, ...
    'control_prediction_error_vs_tracking_error.fig'));
end

function plot_tracking_panel(ax,TP,EX,panel_title)
semilogy(ax,TP.t_set,max(TP.TrackingError_vector,eps),'-', ...
    'LineWidth',1.8,'DisplayName','Test-point online GP'); hold(ax,'on');
semilogy(ax,EX.t_set,max(EX.TrackingError_vector,eps),'--', ...
    'LineWidth',1.8,'DisplayName','Exact dynamics');
grid(ax,'on'); box(ax,'on');
xlabel(ax,'Time (s)'); ylabel(ax,'Tracking error ||\vartheta||');
title(ax,panel_title); legend(ax,'Location','best');
end

function export_clean(fig)
axes_set = findall(fig,'Type','axes');
for ax_i = 1:numel(axes_set)
    disableDefaultInteractivity(axes_set(ax_i));
    set(axes_set(ax_i),'Toolbar',[]);
end
end
