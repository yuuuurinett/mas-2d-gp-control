function OutputFiles = plot_ip_ac_seyboth_simple_diagnostics(do_run)
%PLOT_IP_AC_SEYBOTH_SIMPLE_DIAGNOSTICS Simple 0-30 s AC diagnostics.

if nargin < 1 || isempty(do_run), do_run = true; end
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result','ip_ac_seyboth_m600');
if ~isfolder(result_folder), mkdir(result_folder); end
point_file = fullfile(repo_root,'result','ip_reachable_m600', ...
    'reachable_uniform_M600.mat');
result_name = 'poe_ac_seyboth_R10_T30_simple_diagnostics';
result_file = fullfile(result_folder,[result_name '.mat']);

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_IP_PROJECTION_DIAGNOSTICS','CONTROL_IP_INDUCING_POINT_FILE', ...
    'CONTROL_AC_ITERATION_POLICY','CONTROL_AC_FIXED_ITERATIONS', ...
    'CONTROL_AC_BROADCAST_TRIGGER','CONTROL_AC_SEYBOTH_C0', ...
    'CONTROL_AC_SEYBOTH_C1','CONTROL_AC_SEYBOTH_LAMBDA', ...
    'CONTROL_AC_TRIGGER_DIAGNOSTICS'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment_local( ...
    env_names,old_values));

setenv(env_names{1},'30');
setenv(env_names{2},'600');
setenv(env_names{3},'0');
setenv(env_names{4},point_file);
setenv(env_names{5},'fixed');
setenv(env_names{6},'10');
setenv(env_names{7},'seyboth');
setenv(env_names{8},'0.5');
setenv(env_names{9},'1');
setenv(env_names{10},'5');
setenv(env_names{11},'1');

if do_run || ~isfile(result_file)
    run_simulation_inducing_point('poe_ac',result_folder,result_name, ...
        false,[],[],[],42);
end

d = load(result_file,'t_set','ACOnlineUpdateInterval', ...
    'ACFixedIterations','ac_iteration_count_set', ...
    'ac_inner_trigger_instance_set','ac_trigger_measure_set', ...
    'ac_trigger_threshold_set','projection_absolute_change_rms_set', ...
    'ac_broadcasts_per_agent');
update_idx = find(d.ac_iteration_count_set > 0);
rounds = d.ACFixedIterations;
round_offset = (0:rounds-1).'*(d.ACOnlineUpdateInterval/rounds);
sample_time = round_offset+d.t_set(update_idx);

measure = d.ac_trigger_measure_set(:,1:rounds,update_idx);
network_error = squeeze(max(measure,[],1));
threshold = squeeze(d.ac_trigger_threshold_set(1,1:rounds,update_idx));
event_mask = squeeze(any( ...
    d.ac_inner_trigger_instance_set(:,1:rounds,update_idx),1));
sample_time = sample_time(:);
network_error = network_error(:);
threshold = threshold(:);
event_mask = logical(event_mask(:));

full_file = fullfile(result_folder, ...
    'seyboth_trigger_error_and_inducing_change_0_30s.png');
fig = figure('Visible','off','Color','w','Position',[100 60 1000 720]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
plot(ax1,sample_time,network_error,'Color',[0.85 0.10 0.55], ...
    'LineWidth',0.8,'DisplayName','Maximum local broadcast error');
plot(ax1,sample_time,threshold,'k-','LineWidth',1.0, ...
    'DisplayName','Seyboth threshold');
plot(ax1,sample_time(event_mask),network_error(event_mask),'o', ...
    'Color',[0.85 0.10 0.55],'MarkerFaceColor','none', ...
    'MarkerSize',3,'LineStyle','none','DisplayName','Broadcast event');
xlim(ax1,[0 30]); ylabel(ax1,'RMS error');
title(ax1,sprintf('Seyboth AC trigger, %.1f broadcasts/agent', ...
    d.ac_broadcasts_per_agent));
legend(ax1,'Location','northeast');

ax2 = nexttile(tl); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
p_change = d.projection_absolute_change_rms_set;
valid = isfinite(p_change);
plot(ax2,d.t_set(valid),p_change(valid),'Color',[0.90 0.35 0.05], ...
    'LineWidth',1.2);
xlim(ax2,[0 30]); xlabel(ax2,'Time (s)');
ylabel(ax2,'RMS change in P');
title(ax2,'Change of inducing-point information between updates');
print(fig,full_file,'-dpng','-r180');
close(fig);

zoom_file = fullfile(result_folder, ...
    'seyboth_trigger_error_zoom_10_10p5s.png');
fig = figure('Visible','off','Color','w','Position',[120 100 900 430]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plot(ax,sample_time,network_error,'Color',[0.85 0.10 0.55], ...
    'LineWidth',1.3,'DisplayName','Maximum local broadcast error');
plot(ax,sample_time,threshold,'k-','LineWidth',1.2, ...
    'DisplayName','Seyboth threshold');
plot(ax,sample_time(event_mask),network_error(event_mask),'o', ...
    'Color',[0.85 0.10 0.55],'MarkerFaceColor','none', ...
    'MarkerSize',5,'LineStyle','none','DisplayName','Broadcast event');
xlim(ax,[10 10.5]); xlabel(ax,'Time (s)'); ylabel(ax,'RMS error');
title(ax,'Seyboth AC trigger - enlarged view');
legend(ax,'Location','northeast');
print(fig,zoom_file,'-dpng','-r180');
close(fig);

OutputFiles = string({full_file;zoom_file});
fprintf('Saved simple Seyboth diagnostics:\n  %s\n  %s\n', ...
    full_file,zoom_file);
end

function restore_environment_local(names,values)
for i = 1:numel(names), setenv(names{i},values{i}); end
end
