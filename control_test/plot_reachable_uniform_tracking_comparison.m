function output_file = plot_reachable_uniform_tracking_comparison()
%PLOT_REACHABLE_UNIFORM_TRACKING_COMPARISON Full-trajectory comparison.

repo_root = fileparts(fileparts(mfilename('fullpath')));
base = load(fullfile(repo_root,'result','ip_agent_kappa_sweep', ...
    'poe_M400_R1_eps_0p1_kappa_1_T30.mat'), ...
    't_set','TrackingError_vector');
reach = load(fullfile(repo_root,'result','ip_reachable_uniform', ...
    'poe_M400_reachable_uniform_R1_eps_0p1_kappa_1_T30.mat'), ...
    't_set','TrackingError_vector');

fig = figure('Color','w','Position',[100 100 1050 520]);
ax = axes(fig); hold(ax,'on');
plot(ax,base.t_set,base.TrackingError_vector, ...
    'LineWidth',1.8,'DisplayName','Full-domain uniform');
plot(ax,reach.t_set,reach.TrackingError_vector, ...
    'LineWidth',1.8,'DisplayName','Reachable-box uniform');
xlabel(ax,'Time (s)');
ylabel(ax,'Tracking error norm');
title(ax,'IP-DAC tracking error: inducing-point sampling domain');
xlim(ax,[0,30]);
grid(ax,'on'); box(ax,'on');
legend(ax,'Location','northeast');

result_folder = fullfile(repo_root,'result','ip_reachable_uniform');
output_file = fullfile(result_folder,'tracking_error_comparison.png');
exportgraphics(fig,output_file,'Resolution',240);
savefig(fig,fullfile(result_folder,'tracking_error_comparison.fig'));
end
