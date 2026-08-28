function plot_cen_nbr_tp_exact_every_step(ResultFolder)
% Compare every-step online centralized, neighbor and test-point baselines.

if nargin < 1 || isempty(ResultFolder)
    ProjectRoot = fileparts(fileparts(mfilename('fullpath')));
    ResultFolder = fullfile(ProjectRoot,'result','tp_ip_every_step_baseline');
end

C = load(fullfile(ResultFolder,'centralized_poe_every_step.mat'));
N = load(fullfile(ResultFolder,'neighbor_poe_every_step.mat'));
T = load(fullfile(ResultFolder,'test_point_poe_every_step.mat'));
E = load(fullfile(ResultFolder,'exact_true_dynamics.mat'));

%% Closed-loop tracking comparison
fig = figure('Color','w','Position',[70,70,1100,540]);
semilogy(E.t_set,max(E.TrackingError_vector,eps),'k--','LineWidth',2.0, ...
    'DisplayName','Exact dynamics'); hold on;
semilogy(C.t_set,max(C.TrackingError_vector,eps),'-','LineWidth',1.7, ...
    'DisplayName','Centralized PoE');
semilogy(N.t_set,max(N.TrackingError_vector,eps),'-','LineWidth',1.7, ...
    'DisplayName','Neighbor PoE');
semilogy(T.t_set,max(T.TrackingError_vector,eps),'-','LineWidth',1.7, ...
    'DisplayName','Test-point PoE-DAC');
grid on; box on;
xlabel('Time (s)'); ylabel('Tracking error ||\vartheta||');
title('Tracking error with every-step online learning');
legend('Location','northeast');
disableDefaultInteractivity(gca); set(gca,'Toolbar',[]);
exportgraphics(fig,fullfile(ResultFolder, ...
    'tracking_error_centralized_neighbor_tp_exact.png'),'Resolution',220);
savefig(fig,fullfile(ResultFolder, ...
    'tracking_error_centralized_neighbor_tp_exact.fig'));

%% Prediction-error norms (each method evaluated on its own trajectory)
fig = figure('Color','w','Position',[70,70,1100,540]);
semilogy(C.t_set,max(C.prediction_error_norm_vector,eps),'-', ...
    'LineWidth',1.6,'DisplayName','Centralized PoE'); hold on;
semilogy(N.t_set,max(N.prediction_error_norm_vector,eps),'-', ...
    'LineWidth',1.6,'DisplayName','Neighbor PoE');
semilogy(T.t_set,max(T.prediction_error_norm_vector,eps),'-', ...
    'LineWidth',1.6,'DisplayName','Test-point PoE-DAC');
grid on; box on;
xlabel('Time (s)'); ylabel('Prediction error norm');
title('Unknown-dynamics prediction error (own closed-loop trajectory)');
legend('Location','best');
disableDefaultInteractivity(gca); set(gca,'Toolbar',[]);
exportgraphics(fig,fullfile(ResultFolder, ...
    'prediction_error_centralized_neighbor_tp.png'),'Resolution',220);
savefig(fig,fullfile(ResultFolder, ...
    'prediction_error_centralized_neighbor_tp.fig'));

%% Predicted versus true dynamics, one row per method
methods = {C,N,T};
method_names = {'Centralized PoE','Neighbor PoE','Test-point PoE-DAC'};
y_dim = size(C.f_true_all_set,1);
colors = lines(size(C.f_true_all_set,2));
fig = figure('Color','w','Position',[25,25,1550,1100]);
layout = tiledlayout(fig,numel(methods),y_dim, ...
    'TileSpacing','compact','Padding','compact');
for method_i = 1:numel(methods)
    S = methods{method_i};
    for output_i = 1:y_dim
        ax = nexttile(layout); hold(ax,'on');
        for agent_i = 1:size(S.f_true_all_set,2)
            plot(ax,S.t_set,squeeze(S.f_true_all_set(output_i,agent_i,:)), ...
                '-','Color',colors(agent_i,:),'LineWidth',0.9, ...
                'HandleVisibility','off');
            plot(ax,S.t_set,squeeze(S.f_hat_all_set(output_i,agent_i,:)), ...
                '--','Color',colors(agent_i,:),'LineWidth',0.9, ...
                'DisplayName',sprintf('Agent %d',agent_i));
        end
        grid(ax,'on'); box(ax,'on');
        xlabel(ax,'Time (s)'); ylabel(ax,sprintf('f_%d',output_i));
        title(ax,sprintf('%s: output %d',method_names{method_i},output_i));
        disableDefaultInteractivity(ax); set(ax,'Toolbar',[]);
    end
end
title(layout,'True dynamics (solid) and GP prediction (dashed)');
legend(nexttile(layout,1),'Location','eastoutside');
exportgraphics(fig,fullfile(ResultFolder, ...
    'predicted_vs_true_centralized_neighbor_tp.png'),'Resolution',220);
savefig(fig,fullfile(ResultFolder, ...
    'predicted_vs_true_centralized_neighbor_tp.fig'));
end
