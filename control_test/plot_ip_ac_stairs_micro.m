function output_file = plot_ip_ac_stairs_micro(agent_i, sigma_tag, zoom_start, zoom_span_s)
%PLOT_IP_AC_STAIRS_MICRO  Paper-consistent PETC detector visualization.
% The macro curve uses one trigger ratio per online GP window, avoiding
% misleading vertical connections between successive assignment resets.
% The inset shows the actual squared event condition and distinguishes the
% pre-broadcast detector value from the post-broadcast reset to zero.
%
%   agent_i      : which agent to plot (default 4)
%   sigma_tag    : which sweep result to load, e.g. '0p3' for sigma=0.3
%   zoom_start   : start time (s) of the inset window (default 10.0)
%   zoom_span_s  : inset span in seconds, e.g. 0.1 (1 window) or 0.2
%                  (2 windows) (default 0.2, per the latest request)

if nargin < 1 || isempty(agent_i), agent_i = 4; end
if nargin < 2 || isempty(sigma_tag), sigma_tag = '0p3'; end
if nargin < 3 || isempty(zoom_start), zoom_start = 10.0; end
if nargin < 4 || isempty(zoom_span_s), zoom_span_s = 0.2; end

repo_root = fileparts(fileparts(mfilename('fullpath')));
result_folder = fullfile(repo_root,'result','ip_ac_meng_sigma_sweep');
result_file = fullfile(result_folder, ...
    sprintf('poe_ac_petc_cache_R10_sigma_%s_T30.mat',sigma_tag));

d = load(result_file,'t_set','ac_iteration_count_set', ...
    'ac_inner_trigger_instance_set','ac_trigger_measure_set', ...
    'ac_trigger_threshold_set','ACOnlineUpdateInterval', ...
    'ACFixedIterations');

update_idx = find(d.ac_iteration_count_set > 0);
rounds = d.ACFixedIterations;
round_step = d.ACOnlineUpdateInterval/rounds;
round_offset = (0:rounds-1).'*round_step;
sample_time = round_offset + d.t_set(update_idx);
sample_time = sample_time(:);

measure_matrix = squeeze(d.ac_trigger_measure_set(agent_i,1:rounds,update_idx));
threshold_matrix = squeeze(d.ac_trigger_threshold_set(agent_i,1:rounds,update_idx));
event_matrix = squeeze(d.ac_inner_trigger_instance_set(agent_i,1:rounds,update_idx));

local_measure = double(measure_matrix(:));
trigger_threshold = double(threshold_matrix(:));
triggered = logical(event_matrix(:));

% The Meng et al. detector is e_i^2 > sigma_i^2*eta_hat_i^2.  The saved
% quantities are their square roots, so square them for a literal plot of
% the paper's condition.
error_squared_matrix = double(measure_matrix).^2;
threshold_squared_matrix = double(threshold_matrix).^2;
denominator_floor = max(1e-12,1e-10*max(threshold_squared_matrix(:)));
trigger_ratio_matrix = error_squared_matrix ./ ...
    max(threshold_squared_matrix,denominator_floor);
[~,critical_round] = max(trigger_ratio_matrix,[],1);
critical_index = sub2ind(size(trigger_ratio_matrix),critical_round, ...
    1:size(trigger_ratio_matrix,2));
window_error_squared = error_squared_matrix(critical_index).';
window_threshold_squared = threshold_squared_matrix(critical_index).';
window_triggered = any(logical(event_matrix),1).';
window_time = d.t_set(update_idx).';

% Honest before/after check, printed to console so the number is always
% visible next to the figure (no silent claim of decay either way).
early_mask = sample_time <= 5;
late_mask  = sample_time >= (max(sample_time)-5);
fprintf('Sanity check -- mean |e_i| first 5s: %.4g, last 5s: %.4g (no decay expected)\n', ...
    mean(local_measure(early_mask)), mean(local_measure(late_mask)));

output_file = fullfile(result_folder, ...
    sprintf('ip_ac_stairs_micro_agent%d_sigma_%s.png',agent_i,sigma_tag));

fig = figure('Visible','off','Color','w','Position',[80 60 1150 750]);
ax_main = axes(fig,'Position',[0.09 0.11 0.86 0.80]); hold(ax_main,'on');
grid(ax_main,'on'); box(ax_main,'on');

plot(ax_main,window_time,window_error_squared,'-','Color',[0.15 0.35 0.75], ...
    'LineWidth',1.0,'DisplayName','broadcast error squared (critical round)');
plot(ax_main,window_time,window_threshold_squared,'--', ...
    'Color',[0.10 0.55 0.25],'LineWidth',1.2, ...
    'DisplayName','PETC threshold squared (same round)');
plot(ax_main,window_time(window_triggered), ...
    window_error_squared(window_triggered),'o','Color',[0.78 0.12 0.12], ...
    'MarkerSize',3.2,'LineWidth',0.8,'DisplayName','window contains broadcast');
xlim(ax_main,[0 max(window_time)]);
xlabel(ax_main,'t (s)');
ylabel(ax_main,'Squared detector quantities');
title(ax_main,sprintf(['Agent %d periodic ET-AC: one detector summary per ' ...
    '%.2f s GP update'],agent_i,d.ACOnlineUpdateInterval));
legend(ax_main,'Location','southwest','Interpreter','none');

% --- Inset: stairs-style single/double window zoom with star markers ---
zoom_end = zoom_start + zoom_span_s;
zoom_mask = sample_time >= zoom_start & sample_time < zoom_end;

ax_inset = axes(fig,'Position',[0.58 0.55 0.34 0.32]); hold(ax_inset,'on');
box(ax_inset,'on'); grid(ax_inset,'on');
local_error_squared = local_measure.^2;
threshold_squared = trigger_threshold.^2;
stairs(ax_inset,sample_time(zoom_mask),local_error_squared(zoom_mask),'-', ...
    'Color',[0.15 0.35 0.75],'LineWidth',1.4);
stairs(ax_inset,sample_time(zoom_mask),threshold_squared(zoom_mask),'--', ...
    'Color',[0.10 0.60 0.25],'LineWidth',1.2);
zoom_triggered = zoom_mask & triggered;
plot(ax_inset,sample_time(zoom_triggered),local_error_squared(zoom_triggered),'*', ...
    'Color',[0.85 0.10 0.10],'MarkerSize',9,'LineWidth',1.4);
trigger_times = sample_time(zoom_triggered);
trigger_values = local_error_squared(zoom_triggered);
for event_i = 1:numel(trigger_times)
    plot(ax_inset,[trigger_times(event_i),trigger_times(event_i)], ...
        [trigger_values(event_i),0],':','Color',[0.75 0.15 0.15], ...
        'LineWidth',0.8,'HandleVisibility','off');
end
plot(ax_inset,trigger_times,zeros(size(trigger_times)),'o', ...
    'Color',[0.75 0.15 0.15],'MarkerFaceColor','w','MarkerSize',4, ...
    'LineWidth',0.8);
xlim(ax_inset,[zoom_start, zoom_end]);
title(ax_inset,sprintf('Zoom: [%.2f,%.2f]s (%d rounds)',zoom_start,zoom_end, ...
    nnz(zoom_mask)),'FontSize',9);
xlabel(ax_inset,'t (s)','FontSize',8);
ylabel(ax_inset,'Error squared / threshold squared','FontSize',8);
set(ax_inset,'FontSize',7);

print(fig,output_file,'-dpng','-r200');
savefig(fig,strrep(output_file,'.png','.fig'));
close(fig);

fprintf('Agent %d: %d total triggers over %d windows. Inset [%.2f,%.2f]s: %d/%d rounds triggered.\n', ...
    agent_i,nnz(triggered),numel(update_idx),zoom_start,zoom_end, ...
    nnz(zoom_triggered),nnz(zoom_mask));
fprintf('Saved %s\n',output_file);
end
