function output_file = plot_ip_ac_seyboth_tracking_T50()
%PLOT_IP_AC_SEYBOTH_TRACKING_T50 Simple full and post-transient curves.

repo_root = fileparts(fileparts(mfilename('fullpath')));
result_folder = fullfile(repo_root,'result', ...
    'ip_ac_seyboth_threshold_sweep');
result_file = fullfile(result_folder, ...
    'poe_ac_R10_c0_0p75_c1_1p5_T50.mat');
d = load(result_file,'t_set','TrackingError_vector');

output_file = fullfile(result_folder, ...
    'seyboth_ac_tracking_error_T50.png');
curve_color = [0.85 0.20 0.12];
fig = figure('Visible','off','Color','w','Position',[100 60 950 680]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
plot(ax1,d.t_set,d.TrackingError_vector,'Color',curve_color, ...
    'LineWidth',1.4);
xline(ax1,30,'--','30 s','Color',[0.35 0.35 0.35], ...
    'LabelVerticalAlignment','bottom');
xlim(ax1,[0 50]); ylabel(ax1,'Tracking error');
title(ax1,'IP-AC tracking error over 50 s');

ax2 = nexttile(tl); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
plot(ax2,d.t_set,d.TrackingError_vector,'Color',curve_color, ...
    'LineWidth',1.4);
xline(ax2,30,'--','30 s','Color',[0.35 0.35 0.35], ...
    'LabelVerticalAlignment','bottom');
xlim(ax2,[10 50]); ylim(ax2,[0 0.25]);
xlabel(ax2,'Time (s)'); ylabel(ax2,'Tracking error');
title(ax2,'Post-transient enlarged view');

print(fig,output_file,'-dpng','-r180');
close(fig);

late_mask = d.t_set >= 30;
[late_min,local_min_idx] = min(d.TrackingError_vector(late_mask));
[late_max,local_max_idx] = max(d.TrackingError_vector(late_mask));
late_time = d.t_set(late_mask);
fprintf(['After 30 s: min %.6f at %.2f s; max %.6f at %.2f s; ', ...
    'final %.6f.\n'],late_min,late_time(local_min_idx), ...
    late_max,late_time(local_max_idx),d.TrackingError_vector(end));
fprintf('Saved %s\n',output_file);
end
