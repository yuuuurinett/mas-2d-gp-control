function output_file = plot_ip_ac_diff_curve(agent_i, sigma_tag, t_range)
%PLOT_IP_AC_DIFF_CURVE  Plot g(t) = LHS - RHS as a single curve, per
%window, so the trigger instant is simply where g(t) crosses/reaches 0 -
%avoids the ambiguity of comparing two independently-updating step
%functions (LHS updates every round; RHS only jumps when a neighbor
%broadcasts).
%
%   agent_i  : which agent to plot (default 4)
%   sigma_tag: which sweep result to load, e.g. '0p3' for sigma=0.3
%   t_range  : [t_start t_end] time window to display (default [0 30])

if nargin < 1 || isempty(agent_i), agent_i = 4; end
if nargin < 2 || isempty(sigma_tag), sigma_tag = '0p3'; end
if nargin < 3 || isempty(t_range), t_range = [0 30]; end

repo_root = fileparts(fileparts(mfilename('fullpath')));
result_folder = fullfile(repo_root,'result','ip_ac_meng_sigma_sweep');
result_file = fullfile(result_folder, ...
    sprintf('poe_ac_petc_cache_R10_sigma_%s_T30.mat',sigma_tag));

d = load(result_file,'t_set','ac_iteration_count_set', ...
    'ac_inner_trigger_instance_set','ac_trigger_measure_set', ...
    'ac_trigger_threshold_set','ACOnlineUpdateInterval','ACFixedIterations');

update_idx = find(d.ac_iteration_count_set > 0);
rounds = d.ACFixedIterations;
round_step = d.ACOnlineUpdateInterval/rounds;

measure_all = squeeze(d.ac_trigger_measure_set(agent_i,1:rounds,update_idx));
threshold_all = squeeze(d.ac_trigger_threshold_set(agent_i,1:rounds,update_idx));
event_all = squeeze(d.ac_inner_trigger_instance_set(agent_i,1:rounds,update_idx));

output_file = fullfile(result_folder, ...
    sprintf('ip_ac_diff_curve_agent%d_sigma_%s.png',agent_i,sigma_tag));

fig = figure('Visible','off','Color','w','Position',[80 60 1200 600]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
yline(ax,0,'k-','LineWidth',1.0,'HandleVisibility','off');

n_windows = numel(update_idx);
for w = 1:n_windows
    t0 = d.t_set(update_idx(w));
    tw = t0 + (0:rounds-1)*round_step;
    if tw(end) < t_range(1) || tw(1) > t_range(2), continue; end

    LHS_w = double(measure_all(:,w)).^2;
    RHS_w = double(threshold_all(:,w)).^2;
    g_w = LHS_w - RHS_w;
    trig_w = logical(event_all(:,w));

    stairs(ax,tw,g_w,'-','Color',[0.20 0.20 0.75],'LineWidth',1.0, ...
        'HandleVisibility','off');
    plot(ax,tw(trig_w),g_w(trig_w),'x','Color',[0.85 0.10 0.10], ...
        'MarkerSize',6,'LineWidth',1.4,'HandleVisibility','off');
end

plot(ax,nan,nan,'-','Color',[0.20 0.20 0.75],'LineWidth',1.0, ...
    'DisplayName','g(t) = LHS - RHS');
plot(ax,nan,nan,'x','Color',[0.85 0.10 0.10],'MarkerSize',6,'LineWidth',1.4, ...
    'DisplayName','Broadcast (g \geq 0 detected)');
yline(ax,0,'k-','LineWidth',1.0,'DisplayName','g=0 (trigger line)');

xlabel(ax,'t (s)');
ylabel(ax,'g(t) = LHS - RHS');
title(ax,sprintf('Agent %d: single difference curve, trigger = crossing 0 (\\sigma=%s)', ...
    agent_i,strrep(sigma_tag,'p','.')));
legend(ax,'Location','best');
xlim(ax,t_range);

print(fig,output_file,'-dpng','-r200');
savefig(fig,strrep(output_file,'.png','.fig'));
close(fig);
fprintf('Saved %s\n',output_file);
end