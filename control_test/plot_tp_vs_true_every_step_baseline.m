function plot_tp_vs_true_every_step_baseline(ResultFolder)
% Test-point baseline against exact-model control and true dynamics.

if nargin < 1 || isempty(ResultFolder)
    ProjectRoot = fileparts(fileparts(mfilename('fullpath')));
    ResultFolder = fullfile(ProjectRoot,'result','tp_ip_every_step_baseline');
end

TP = load(fullfile(ResultFolder,'test_point_poe_every_step.mat'));
EX = load(fullfile(ResultFolder,'exact_true_dynamics.mat'));

%% Tracking: test-point GP controller versus exact-dynamics controller
fig = figure('Color','w','Position',[80,80,1050,520]);
semilogy(TP.t_set,max(TP.TrackingError_vector,eps),'-', ...
    'LineWidth',1.8,'DisplayName','Test-point GP');
hold on;
semilogy(EX.t_set,max(EX.TrackingError_vector,eps),'--', ...
    'LineWidth',1.8,'DisplayName','Exact true dynamics');
grid on; box on;
xlabel('Time (s)'); ylabel('Tracking error ||\vartheta||');
title('Tracking error: test-point GP versus exact dynamics');
legend('Location','northeast');
disableDefaultInteractivity(gca); set(gca,'Toolbar',[]);
exportgraphics(fig,fullfile(ResultFolder, ...
    'tracking_error_test_point_vs_exact.png'),'Resolution',220);
savefig(fig,fullfile(ResultFolder,'tracking_error_test_point_vs_exact.fig'));

%% On the TP trajectory: TP prediction versus true unknown dynamics
AgentQuantity = size(TP.f_true_all_set,2);
y_dim = size(TP.f_true_all_set,1);
colors = lines(AgentQuantity);
fig = figure('Color','w','Position',[50,50,1500,720]);
layout = tiledlayout(fig,1,y_dim,'TileSpacing','compact','Padding','compact');
for output_i = 1:y_dim
    ax = nexttile(layout);
    hold(ax,'on');
    for agent_i = 1:AgentQuantity
        plot(ax,TP.t_set,squeeze(TP.f_true_all_set(output_i,agent_i,:)), ...
            '-','Color',colors(agent_i,:),'LineWidth',1.1, ...
            'DisplayName',sprintf('Agent %d true',agent_i));
        plot(ax,TP.t_set,squeeze(TP.f_hat_all_set(output_i,agent_i,:)), ...
            '--','Color',colors(agent_i,:),'LineWidth',1.1, ...
            'DisplayName',sprintf('Agent %d TP prediction',agent_i));
    end
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'Time (s)'); ylabel(ax,sprintf('f_%d',output_i));
    title(ax,sprintf('Unknown dynamics output %d',output_i));
    disableDefaultInteractivity(ax); set(ax,'Toolbar',[]);
end
title(layout,'Test-point prediction versus true unknown dynamics (solid: true, dashed: TP)');
legend(nexttile(layout,1),'Location','eastoutside');
exportgraphics(fig,fullfile(ResultFolder, ...
    'test_point_prediction_vs_true_dynamics.png'),'Resolution',220);
savefig(fig,fullfile(ResultFolder, ...
    'test_point_prediction_vs_true_dynamics.fig'));
end
