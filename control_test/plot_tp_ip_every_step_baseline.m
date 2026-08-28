function plot_tp_ip_every_step_baseline(ResultFolder)
% Plot the fair every-step-learning TP/IP baseline comparison.

if nargin < 1 || isempty(ResultFolder)
    ProjectRoot = fileparts(fileparts(mfilename('fullpath')));
    ResultFolder = fullfile(ProjectRoot,'result','tp_ip_every_step_baseline');
end

TP = load(fullfile(ResultFolder,'test_point_poe_every_step.mat'));
IP = load(fullfile(ResultFolder,'inducing_point_poe_every_step.mat'));

%% Tracking error
fig = figure('Color','w','Position',[80,80,1050,520]);
semilogy(TP.t_set,max(TP.TrackingError_vector,eps),'-', ...
    'LineWidth',1.8,'DisplayName','Test-point baseline (full communication)');
hold on;
semilogy(IP.t_set,max(IP.TrackingError_vector,eps),'-', ...
    'LineWidth',1.8,'DisplayName','Inducing-point ET (equation 17)');
grid on; box on;
xlabel('Time (s)');
ylabel('Tracking error ||\vartheta||');
title('Tracking error: every-step online learning');
legend('Location','northeast');
disableDefaultInteractivity(gca);
set(gca,'Toolbar',[]);
exportgraphics(fig,fullfile(ResultFolder, ...
    'tracking_error_tp_vs_ip_every_step.png'),'Resolution',220);
savefig(fig,fullfile(ResultFolder,'tracking_error_tp_vs_ip_every_step.fig'));

%% Predicted versus true unknown dynamics
AgentQuantity = size(TP.f_true_all_set,2);
y_dim = size(TP.f_true_all_set,1);
colors = lines(AgentQuantity);

fig = figure('Color','w','Position',[30,30,1500,900]);
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

plot_prediction_panel(nexttile(layout),TP.t_set,TP.f_true_all_set, ...
    TP.f_hat_all_set,1,colors,'Test-point: output 1');
plot_prediction_panel(nexttile(layout),TP.t_set,TP.f_true_all_set, ...
    TP.f_hat_all_set,min(2,y_dim),colors,'Test-point: output 2');
plot_prediction_panel(nexttile(layout),IP.t_set,IP.f_true_all_set, ...
    IP.f_hat_all_set,1,colors,'Inducing-point ET: output 1');
plot_prediction_panel(nexttile(layout),IP.t_set,IP.f_true_all_set, ...
    IP.f_hat_all_set,min(2,y_dim),colors,'Inducing-point ET: output 2');

title(layout,['Predicted and true unknown dynamics ', ...
    '(solid: true, dashed: prediction)']);
legend(nexttile(layout,1),'Location','eastoutside');
all_axes = findall(fig,'Type','axes');
for ax_i = 1:numel(all_axes)
    disableDefaultInteractivity(all_axes(ax_i));
    set(all_axes(ax_i),'Toolbar',[]);
end
exportgraphics(fig,fullfile(ResultFolder, ...
    'predicted_vs_true_dynamics_tp_vs_ip_every_step.png'),'Resolution',220);
savefig(fig,fullfile(ResultFolder, ...
    'predicted_vs_true_dynamics_tp_vs_ip_every_step.fig'));

%% Prediction-error norm, useful for reading the dynamics comparison
fig = figure('Color','w','Position',[80,80,1050,520]);
semilogy(TP.t_set,max(TP.prediction_error_norm_vector,eps),'-', ...
    'LineWidth',1.6,'DisplayName','Test-point baseline');
hold on;
semilogy(IP.t_set,max(IP.prediction_error_norm_vector,eps),'-', ...
    'LineWidth',1.6,'DisplayName','Inducing-point ET');
grid on; box on;
xlabel('Time (s)');
ylabel('Prediction error norm');
title('Unknown-dynamics prediction error');
legend('Location','best');
disableDefaultInteractivity(gca);
set(gca,'Toolbar',[]);
exportgraphics(fig,fullfile(ResultFolder, ...
    'prediction_error_tp_vs_ip_every_step.png'),'Resolution',220);
savefig(fig,fullfile(ResultFolder,'prediction_error_tp_vs_ip_every_step.fig'));
end

function plot_prediction_panel(ax,t,f_true,f_hat,output_idx,colors,panel_title)
hold(ax,'on');
AgentQuantity = size(f_true,2);
for agent_i = 1:AgentQuantity
    plot(ax,t,squeeze(f_true(output_idx,agent_i,:)),'-', ...
        'Color',colors(agent_i,:),'LineWidth',1.0, ...
        'DisplayName',sprintf('Agent %d true',agent_i));
    plot(ax,t,squeeze(f_hat(output_idx,agent_i,:)),'--', ...
        'Color',colors(agent_i,:),'LineWidth',1.0, ...
        'DisplayName',sprintf('Agent %d prediction',agent_i));
end
grid(ax,'on'); box(ax,'on');
xlabel(ax,'Time (s)'); ylabel(ax,sprintf('f_%d',output_idx));
title(ax,panel_title);
end
