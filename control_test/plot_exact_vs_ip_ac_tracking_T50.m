function output_file = plot_exact_vs_ip_ac_tracking_T50()
%PLOT_EXACT_VS_IP_AC_TRACKING_T50 Compare like-for-like 50 s tracking errors.

repo_root = fileparts(fileparts(mfilename('fullpath')));
result_folder = fullfile(repo_root,'result', ...
    'ip_ac_seyboth_threshold_sweep');
ip_file = fullfile(result_folder, ...
    'poe_ac_R10_c0_0p75_c1_1p5_T50.mat');
exact_file = fullfile(result_folder,'exact_T50.mat');

ip = load(ip_file,'t_set','TrackingError_vector');
ex = load(exact_file,'t_set','TrackingError_vector');

output_file = fullfile(result_folder, ...
    'tracking_error_exact_vs_ip_ac_T50.png');
fig = figure('Visible','off','Color','w','Position',[100 80 980 560]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plot(ax,ex.t_set,ex.TrackingError_vector,'Color',[0.15 0.15 0.15], ...
    'LineWidth',1.7,'DisplayName','Exact dynamics');
plot(ax,ip.t_set,ip.TrackingError_vector,'Color',[0.85 0.20 0.12], ...
    'LineWidth',1.45,'DisplayName','IP-AC (R=10)');
xlim(ax,[0 50]);
xlabel(ax,'Time (s)'); ylabel(ax,'Tracking error');
title(ax,'Tracking error: Exact dynamics vs. IP-AC');
legend(ax,'Location','northeast');
print(fig,output_file,'-dpng','-r200');
savefig(fig,strrep(output_file,'.png','.fig'));
close(fig);

exact_late = ex.TrackingError_vector(ex.t_set >= 30);
ip_late = ip.TrackingError_vector(ip.t_set >= 30);
fprintf('Exact: mean(full)=%.6f, mean(30-50)=%.6f, final=%.6f\n', ...
    mean(ex.TrackingError_vector),mean(exact_late), ...
    ex.TrackingError_vector(end));
fprintf('IP-AC: mean(full)=%.6f, mean(30-50)=%.6f, final=%.6f\n', ...
    mean(ip.TrackingError_vector),mean(ip_late), ...
    ip.TrackingError_vector(end));
fprintf('Saved %s\n',output_file);
end
