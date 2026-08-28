function output_file = plot_ip_ac_diff_curve(agent_i, sigma_tag, t_range)
%PLOT_IP_AC_DIFF_CURVE_CLEAN  g(t)=LHS-RHS as a thin linearly-interpolated
%zigzag (matching the Umlauft/Meng/Xu reference-figure style) instead of
%a stairs plot, so 300 windows over 30s render as a readable sawtooth
%line rather than a dense filled band.
%
%   agent_i  : which agent to plot (default 4)
%   sigma_tag: which sweep result to load, e.g. '0p3' for sigma=0.3
%   t_range  : [t_start t_end] time window to display (default [0 30])

if nargin < 1 || isempty(agent_i), agent_i = 4; end
if nargin < 2 || isempty(sigma_tag), sigma_tag = '0p3'; end
if nargin < 3 || isempty(t_range), t_range = [0 30]; end

repo_root = fileparts(fileparts(mfilename('fullpath')));
result_folder = fullfile(repo_root,'result', ...
    'ip_ac_meng_sigma_sweep_m400_online300_per_agent');
result_file = fullfile(result_folder, ...
    sprintf('poe_ac_petc_R10_M400_online300_per_agent_sigma_%s_T30.mat',sigma_tag));

d = load(result_file,'t_set','ac_iteration_count_set', ...
    'ac_inner_trigger_instance_set','ac_trigger_measure_set', ...
    'ac_trigger_threshold_set','ACOnlineUpdateInterval','ACFixedIterations');

update_idx = find(d.ac_iteration_count_set > 0);
rounds = d.ACFixedIterations;
round_step = d.ACOnlineUpdateInterval/rounds;

measure_all = squeeze(d.ac_trigger_measure_set(agent_i,1:rounds,update_idx));
threshold_all = squeeze(d.ac_trigger_threshold_set(agent_i,1:rounds,update_idx));
event_all = squeeze(d.ac_inner_trigger_instance_set(agent_i,1:rounds,update_idx));

% Build one long vector across ALL windows, inserting a NaN gap between
% windows so MATLAB's plot() breaks the line there instead of drawing a
% spurious connector -- but WITHIN a window, points are linearly
% interpolated (a thin zigzag), not stepped.
n_windows = numel(update_idx);
t_all = [];
g_all = [];
trig_t = [];
trig_g = [];
for w = 1:n_windows
    t0 = d.t_set(update_idx(w));
    tw = t0 + (0:rounds-1)*round_step;
    if tw(end) < t_range(1) || tw(1) > t_range(2), continue; end

    % The simulation already stores LHS and RHS as normalized Frobenius
    % norms. Do not square them again: the plotted quantity is exactly the
    % detector used by the code, g=LHS-RHS.
    LHS_w = double(measure_all(:,w));
    RHS_w = double(threshold_all(:,w));
    g_w = LHS_w - RHS_w;
    trig_w = logical(event_all(:,w));

    t_all = [t_all, tw, NaN]; %#ok<AGROW>
    g_all = [g_all, g_w(:).', NaN]; %#ok<AGROW>
    trig_t = [trig_t, tw(trig_w)]; %#ok<AGROW>
    trig_g = [trig_g, g_w(trig_w).']; %#ok<AGROW>
end

output_file = fullfile(result_folder, ...
    sprintf('ip_ac_diff_curve_clean_agent%d_sigma_%s.png',agent_i,sigma_tag));

fig = figure('Visible','off','Color','w','Position',[40 40 2400 700]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');

yline(ax,0,'k-','LineWidth',1.0,'DisplayName','g=0 (trigger line)');
plot(ax,t_all,g_all,'-','Color',[0.15 0.30 0.75],'LineWidth',0.35, ...
    'DisplayName','g(t) = LHS - RHS');
plot(ax,trig_t,trig_g,'x','Color',[0.85 0.10 0.10],'MarkerSize',4, ...
    'LineWidth',1.0,'DisplayName','Broadcast (g \geq 0 detected)');

xlabel(ax,'t (s)');
ylabel(ax,'g(t) = LHS - RHS');
title(ax,sprintf('Agent %d: difference curve (linear zigzag, \\sigma=%s)', ...
    agent_i,strrep(sigma_tag,'p','.')));
legend(ax,'Location','best');
xlim(ax,t_range);

print(fig,output_file,'-dpng','-r220');
savefig(fig,strrep(output_file,'.png','.fig'));
close(fig);
fprintf('Saved %s\n',output_file);
end
