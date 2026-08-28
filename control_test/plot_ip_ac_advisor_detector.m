function output_file = plot_ip_ac_advisor_detector( ...
    result_file,agent_i,zoom_start,zoom_duration)
%PLOT_IP_AC_ADVISOR_DETECTOR Exact LHS/RHS trigger history for IP-AC.
%
% The physical GP update index k and the inner AC round ell are flattened
% onto one time axis:
%       t_check = t_k + (ell-1)*(0.1/R).
% At every check the saved pre-broadcast LHS and RHS are plotted. If
% LHS>=RHS triggers a broadcast, XiHat_i is replaced by Xi_i and the local
% broadcast error (LHS) is exactly zero immediately after that check.

if nargin < 2 || isempty(agent_i), agent_i = 4; end
if nargin < 3 || isempty(zoom_start), zoom_start = 10; end
if nargin < 4 || isempty(zoom_duration), zoom_duration = 0.2; end

d = load(result_file,'t_set','ac_iteration_count_set', ...
    'ac_inner_trigger_instance_set','ac_trigger_measure_set', ...
    'ac_trigger_threshold_set','ac_tracking_error_history_set', ...
    'ACFixedIterations','ACPeriodicSigma','NumInducingPoints', ...
    'ACOnlineUpdateInterval');

required_fields = {'ac_trigger_measure_set','ac_trigger_threshold_set', ...
    'ac_inner_trigger_instance_set'};
for field_i = 1:numel(required_fields)
    if ~isfield(d,required_fields{field_i}) || isempty(d.(required_fields{field_i}))
        error(['Result file does not contain AC trigger diagnostics. ' ...
            'Rerun with CONTROL_AC_TRIGGER_DIAGNOSTICS=1.']);
    end
end

rounds = d.ACFixedIterations;
update_idx = find(d.ac_iteration_count_set > 0);
if isempty(update_idx)
    error('No AC update windows were stored in %s.',result_file);
end
window_time = reshape(d.t_set(update_idx),1,[]);
window_count = numel(update_idx);

lhs_matrix = double(squeeze( ...
    d.ac_trigger_measure_set(agent_i,1:rounds,update_idx)));
rhs_matrix = double(squeeze( ...
    d.ac_trigger_threshold_set(agent_i,1:rounds,update_idx)));
event_matrix = logical(squeeze( ...
    d.ac_inner_trigger_instance_set(agent_i,1:rounds,update_idx)));
if isscalar(window_time)
    lhs_matrix = lhs_matrix(:);
    rhs_matrix = rhs_matrix(:);
    event_matrix = event_matrix(:);
end

%% Flatten every real AC check onto the 0--30 s time axis.
round_step = d.ACOnlineUpdateInterval/rounds;
round_offset = (0:rounds-1).'*round_step;
check_time_matrix = round_offset+window_time;
check_time = check_time_matrix(:).';
lhs_pre = lhs_matrix(:).';
rhs_check = rhs_matrix(:).';
triggered = event_matrix(:).';

% Descriptive 0.1 s-window summaries for the long-time overview. These
% means are never used by the detector; the detector still evaluates every
% original (k,ell) pair stored above.
lhs_round_mean = mean(lhs_matrix,1,'omitnan');
rhs_round_mean = mean(rhs_matrix,1,'omitnan');
events_per_window = sum(event_matrix,1);
window_has_event = events_per_window > 0;

% Show the exact post-broadcast reset without pretending that it is a new
% periodic detector check. The short offset is only a plotting phase inside
% the same communication round.
post_offset = 0.22*round_step;
phase_time = reshape([check_time;check_time+post_offset],1,[]);
lhs_post = lhs_pre;
lhs_post(triggered) = 0;
phase_lhs = reshape([lhs_pre;lhs_post],1,[]);
phase_rhs = reshape([rhs_check;rhs_check],1,[]);

zoom_end = zoom_start+zoom_duration;
zoom_phase_mask = phase_time >= zoom_start & phase_time <= zoom_end;
zoom_trigger_mask = triggered & check_time >= zoom_start & ...
    check_time <= zoom_end;

%% AC-error reduction for the first complete 0.1 s window in the zoom.
[~,zoom_window] = min(abs(window_time-zoom_start));
zoom_window_time = window_time(zoom_window);
zoom_error_start = NaN;
zoom_error_end = NaN;
if isfield(d,'ac_tracking_error_history_set') && ...
        numel(d.ac_tracking_error_history_set) >= update_idx(zoom_window)
    error_history = d.ac_tracking_error_history_set{update_idx(zoom_window)};
    if ~isempty(error_history)
        zoom_error_start = error_history(1);
        zoom_error_end = error_history(end);
    end
end

result_folder = fileparts(result_file);
output_file = fullfile(result_folder,sprintf( ...
    'ip_ac_lhs_rhs_round_mean_agent%d_zoom_%0.1fs.png', ...
    agent_i,zoom_start));
output_pdf = strrep(output_file,'.png','.pdf');
exact_output_file = fullfile(result_folder,sprintf( ...
    'ip_ac_lhs_rhs_flattened_checks_agent%d_zoom_%0.1fs.png', ...
    agent_i,zoom_start));
exact_output_pdf = strrep(exact_output_file,'.png','.pdf');

lhs_color = [0.24 0.18 0.48];
rhs_color = [0.78 0.43 0.10];
event_color = [0.72 0.08 0.12];

% Three times the previous horizontal canvas: 3000 periodic checks need
% enough pixels to remain visually separated in the 0--30 s panel.
fig = figure('Visible','off','Color','w','Position',[40 60 3800 1050]);
layout = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
layout.OuterPosition = [0.035 0.095 0.93 0.84];

%% Full trajectory: descriptive mean across the ten AC rounds.
% This panel is deliberately a readable long-time summary. It does not
% replace the R trigger checks and is not used to decide a broadcast. The
% exact per-round detector and reset sequence remains in the zoom panel.
ax_main = nexttile(layout); hold(ax_main,'on'); box(ax_main,'on');
grid(ax_main,'on');
plot(ax_main,window_time,lhs_round_mean,'-', ...
    'Color',lhs_color,'LineWidth',1.25, ...
    'DisplayName',sprintf('Mean LHS over R=%d rounds',rounds));
plot(ax_main,window_time,rhs_round_mean,'--', ...
    'Color',rhs_color,'LineWidth',1.25, ...
    'DisplayName',sprintf('Mean RHS over R=%d rounds',rounds));
plot(ax_main,window_time(window_has_event), ...
    lhs_round_mean(window_has_event),'*', ...
    'LineStyle','none','Color',event_color, ...
    'MarkerSize',3.2,'LineWidth',0.5, ...
    'DisplayName','Window containing actual broadcast(s)');
xlim(ax_main,[window_time(1), ...
    window_time(end)+d.ACOnlineUpdateInterval]);
xticks(ax_main,0:2:ceil(window_time(end)));
xlabel(ax_main,'Time (s)');
ylabel(ax_main,'Detector value');
title(ax_main,sprintf(['0--30 s round mean (display only); detector used ' ...
    'all %d checks and produced %d broadcasts'], ...
    numel(check_time),nnz(triggered)), ...
    'FontWeight','normal');
legend(ax_main,'Location','best','Box','off');

%% Fixed physical-time zoom showing the actual discrete causal sequence.
ax_zoom = nexttile(layout); hold(ax_zoom,'on'); box(ax_zoom,'on');
grid(ax_zoom,'on');
plot(ax_zoom,phase_time(zoom_phase_mask),phase_lhs(zoom_phase_mask),'-', ...
    'Color',lhs_color,'LineWidth',1.45, ...
    'DisplayName','LHS (including reset to zero)');
plot(ax_zoom,phase_time(zoom_phase_mask),phase_rhs(zoom_phase_mask),'--', ...
    'Color',rhs_color,'LineWidth',1.25, ...
    'DisplayName','RHS');
plot(ax_zoom,check_time(zoom_trigger_mask),lhs_pre(zoom_trigger_mask),'*', ...
    'Color',event_color,'MarkerSize',8,'LineWidth',1.2, ...
    'DisplayName','Broadcast: LHS >= RHS');
xlim(ax_zoom,[zoom_start,zoom_end]);
xlabel(ax_zoom,'Time (s)');
ylabel(ax_zoom,'Detector value');
if isfinite(zoom_error_start) && isfinite(zoom_error_end)
    zoom_title = sprintf(['Zoom %.2f--%.2f s; window %.2f s AC error ' ...
        '%.3g -> %.3g'],zoom_start,zoom_end,zoom_window_time, ...
        zoom_error_start,zoom_error_end);
else
    zoom_title = sprintf('Zoom %.2f--%.2f s',zoom_start,zoom_end);
end
title(ax_zoom,zoom_title,'FontWeight','normal');
legend(ax_zoom,'Location','best','Box','off');

title(layout,sprintf(['Agent %d periodic ET-AC detector ' ...
    '(M=%d, R=%d, \\sigma=%.1f)'],agent_i,d.NumInducingPoints, ...
    rounds,d.ACPeriodicSigma),'FontWeight','bold');
figure_note = annotation(fig,'textbox',[0.04 0.012 0.92 0.060], ...
    'String',['Note: the 0--30 s panel is the mean of the ten AC ' ...
    'rounds for visualization only. The detector and broadcast count ' ...
    'use every unaveraged round. A star here marks a 0.1 s window ' ...
    'containing at least one broadcast; the zoom shows exact events.'], ...
    'Interpreter','none','EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',8.5,'Color',[0.25 0.25 0.25]);
exportgraphics(fig,output_file,'Resolution',220);
exportgraphics(fig,output_pdf,'ContentType','vector');
savefig(fig,strrep(output_file,'.png','.fig'));

%% Companion figure: all checks in their actual within-window order.
% There is one LHS curve and one RHS curve.  The ten AC rounds belonging
% to a 0.1 s physical window occupy ten consecutive points on this axis.
cla(ax_main); hold(ax_main,'on'); box(ax_main,'on'); grid(ax_main,'on');
plot(ax_main,check_time,lhs_pre,'-', ...
    'Color',lhs_color,'LineWidth',0.60, ...
    'DisplayName','LHS: 3000 ordered checks');
plot(ax_main,check_time,rhs_check,'--', ...
    'Color',rhs_color,'LineWidth',0.60, ...
    'DisplayName','RHS: 3000 ordered checks');
plot(ax_main,check_time(triggered),lhs_pre(triggered),'*', ...
    'LineStyle','none','Color',event_color, ...
    'MarkerSize',2.5,'LineWidth',0.4, ...
    'DisplayName','Every actual broadcast');
xlim(ax_main,[window_time(1), ...
    window_time(end)+d.ACOnlineUpdateInterval]);
xticks(ax_main,0:2:ceil(window_time(end)));
xlabel(ax_main,'Time (s)');
ylabel(ax_main,'Detector value');
title(ax_main,sprintf(['0--30 s ordered detector history: %d physical ' ...
    'windows x %d rounds = %d checks; %d broadcasts'], ...
    window_count,rounds,numel(check_time),nnz(triggered)), ...
    'FontWeight','normal');
legend(ax_main,'Location','best','Box','off');
figure_note.String = ['Note: no trigger data are averaged in the ' ...
    '0--30 s panel. Each 0.1 s physical window contributes ten ' ...
    'consecutive AC checks to the single LHS and RHS curves; every star ' ...
    'is one broadcast produced by the original LHS >= RHS test.'];
exportgraphics(fig,exact_output_file,'Resolution',220);
exportgraphics(fig,exact_output_pdf,'ContentType','vector');
savefig(fig,strrep(exact_output_file,'.png','.fig'));
close(fig);

fprintf(['Saved mean overview: %s\nSaved exact-check version: %s\n' ...
    'Agent %d: %d broadcasts over %d actual checks; zoom contains %d ' ...
    'broadcasts.\nVector copies: %s and %s\n'], ...
    output_file,exact_output_file,agent_i,nnz(triggered),numel(triggered), ...
    nnz(zoom_trigger_mask),output_pdf,exact_output_pdf);
end
