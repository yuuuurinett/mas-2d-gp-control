function Summary = run_ip_ac_reference_gain_tradeoff(do_run)
%RUN_IP_AC_REFERENCE_GAIN_TRADEOFF M600, R10, T10 gain/reference test.
% The chosen gain should keep the ten-round error near/below one while
% reducing the detector overshoot at LHS >= RHS broadcasts.

if nargin < 1 || isempty(do_run), do_run = true; end
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result', ...
    'ip_ac_reference_gain_tradeoff_m600_T10');
if ~isfolder(result_folder), mkdir(result_folder); end

gains = [0.05;0.075;0.10;0.125;0.14;0.15;0.16;0.175;0.20];
case_count = numel(gains);
files = strings(case_count,1);

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_IP_INDUCING_POINT_FILE','CONTROL_AC_FIXED_ITERATIONS', ...
    'CONTROL_AC_PERIODIC_SIGMA','CONTROL_AC_TRIGGER_DIAGNOSTICS', ...
    'CONTROL_ONLINE_AGENT_POLICY','CONTROL_AC_CONSENSUS_STEP'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment_local(env_names,old_values)); %#ok<NASGU>

setenv(env_names{1},'10');
setenv(env_names{2},'600');
setenv(env_names{3},'');
setenv(env_names{4},'10');
setenv(env_names{5},'0.3');
setenv(env_names{6},'1');
setenv(env_names{7},'all_agents');

for case_i = 1:case_count
    setenv(env_names{8},num2str(gains(case_i),'%.12g'));
    tag = strrep(num2str(gains(case_i),'%.3f'),'.','p');
    name = sprintf('poe_ac_M600_R10_gain_%s_T10_seed42',tag);
    files(case_i) = fullfile(result_folder,[name '.mat']);
    if do_run || ~isfile(files(case_i))
        run_simulation_inducing_point('poe_ac',result_folder,name, ...
            false,[],[],[],42);
    end
end

StartErrorMean = nan(case_count,1);
TerminalErrorMean = nan(case_count,1);
TerminalErrorMax = nan(case_count,1);
TerminalBelowOnePercent = nan(case_count,1);
TriggerGapMean = nan(case_count,1);
TriggerGapMedian = nan(case_count,1);
TriggerRelativeGapMedian = nan(case_count,1);
BroadcastsPerAgent = nan(case_count,1);
selected_histories = cell(case_count,1);
selected_time = 5;

fig = figure('Visible','off','Color','w','Position',[80 80 1350 760]);
layout = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax_time = nexttile(layout); hold(ax_time,'on'); box(ax_time,'on'); grid(ax_time,'on');
ax_round = nexttile(layout); hold(ax_round,'on'); box(ax_round,'on'); grid(ax_round,'on');
colors = lines(case_count);

for case_i = 1:case_count
    d = load(files(case_i),'t_set','ac_iteration_count_set', ...
        'ac_tracking_error_history_set','ac_trigger_measure_set', ...
        'ac_trigger_threshold_set','ac_inner_trigger_instance_set', ...
        'ac_event_count_per_agent');
    update_idx = find(d.ac_iteration_count_set > 0);
    start_error = nan(1,numel(update_idx));
    terminal_error = nan(1,numel(update_idx));
    for window_i = 1:numel(update_idx)
        history = d.ac_tracking_error_history_set{update_idx(window_i)};
        start_error(window_i) = history(1);
        terminal_error(window_i) = history(end);
    end

    StartErrorMean(case_i) = mean(start_error,'omitnan');
    TerminalErrorMean(case_i) = mean(terminal_error,'omitnan');
    TerminalErrorMax(case_i) = max(terminal_error,[],'omitnan');
    TerminalBelowOnePercent(case_i) = 100*mean(terminal_error <= 1);
    BroadcastsPerAgent(case_i) = mean(d.ac_event_count_per_agent);

    lhs = double(d.ac_trigger_measure_set(:,:,update_idx));
    rhs = double(d.ac_trigger_threshold_set(:,:,update_idx));
    events = logical(d.ac_inner_trigger_instance_set(:,:,update_idx));
    gap = lhs(events)-rhs(events);
    rhs_event = rhs(events);
    positive_rhs = rhs_event > 1e-8;
    TriggerGapMean(case_i) = mean(gap,'omitnan');
    TriggerGapMedian(case_i) = median(gap,'omitnan');
    TriggerRelativeGapMedian(case_i) = median( ...
        gap(positive_rhs)./rhs_event(positive_rhs),'omitnan');

    window_time = d.t_set(update_idx);
    [~,selected_i] = min(abs(window_time-selected_time));
    selected_histories{case_i} = ...
        d.ac_tracking_error_history_set{update_idx(selected_i)};

    plot(ax_time,window_time,terminal_error,'-', ...
        'Color',colors(case_i,:),'LineWidth',1.1, ...
        'DisplayName',sprintf('h_{AC}=%.3f',gains(case_i)));
    plot(ax_round,0:10,selected_histories{case_i},'-o', ...
        'Color',colors(case_i,:),'MarkerSize',3.5,'LineWidth',1.25, ...
        'DisplayName',sprintf('h_{AC}=%.3f: %.2f -> %.2f', ...
        gains(case_i),selected_histories{case_i}(1), ...
        selected_histories{case_i}(end)));
end

yline(ax_time,1,'k--','Target = 1','LineWidth',1.1, ...
    'HandleVisibility','off');
xlabel(ax_time,'Physical time (s)');
ylabel(ax_time,'RMS error after R=10');
title(ax_time,'Reference-tracking error remaining after ten AC rounds', ...
    'FontWeight','normal');
legend(ax_time,'Location','best','Box','off');

yline(ax_round,1,'k--','Target = 1','LineWidth',1.1, ...
    'HandleVisibility','off');
xlim(ax_round,[0 10]); xticks(ax_round,0:10);
xlabel(ax_round,'AC round');
ylabel(ax_round,'RMS error to current average P');
title(ax_round,sprintf(['Selected %.1f--%.1f s window: ' ...
    'what ten AC rounds accomplish'],selected_time,selected_time+0.1), ...
    'FontWeight','normal');
legend(ax_round,'Location','best','Box','off');
title(layout,'M600 IP-AC consensus-gain trade-off (R=10, sigma=0.3)', ...
    'FontWeight','bold');

plot_file = fullfile(result_folder,'reference_error_gain_tradeoff.png');
exportgraphics(fig,plot_file,'Resolution',230);
exportgraphics(fig,strrep(plot_file,'.png','.pdf'),'ContentType','vector');
savefig(fig,strrep(plot_file,'.png','.fig'));
close(fig);

% A clean single-window view matching the original diagnostic figure.
recommended_i = find(abs(gains-0.10) < 1e-12,1);
recommended_history = selected_histories{recommended_i};
simple_fig = figure('Visible','off','Color','w', ...
    'Position',[120 120 650 450]);
simple_ax = axes(simple_fig); hold(simple_ax,'on');
box(simple_ax,'on'); grid(simple_ax,'on');
plot(simple_ax,0:10,recommended_history,'-o', ...
    'Color',[0.15 0.35 0.65],'MarkerFaceColor','w', ...
    'MarkerSize',4.5,'LineWidth',1.6);
yline(simple_ax,1,'--','target = 1','Color',[0.25 0.25 0.25], ...
    'LineWidth',1.0,'LabelHorizontalAlignment','left');
xlim(simple_ax,[0 10]); xticks(simple_ax,0:10);
ylim(simple_ax,[0,max(4,1.08*recommended_history(1))]);
xlabel(simple_ax,'AC round');
ylabel(simple_ax,'RMS error to current average P');
title(simple_ax,sprintf(['Selected window %.1f--%.1f s: ' ...
    'what ten AC rounds accomplish'],selected_time,selected_time+0.1), ...
    'FontWeight','normal');
text(simple_ax,9.75,0.88*recommended_history(1), ...
    sprintf('%.2f  -->  %.2f\n(%.1f%% reduction)', ...
    recommended_history(1),recommended_history(end), ...
    100*(1-recommended_history(end)/recommended_history(1))), ...
    'HorizontalAlignment','right','VerticalAlignment','top', ...
    'Color',[0.75 0.12 0.12]);
simple_file = fullfile(result_folder, ...
    'reference_error_selected_window_h0p10.png');
exportgraphics(simple_fig,simple_file,'Resolution',230);
exportgraphics(simple_fig,strrep(simple_file,'.png','.pdf'), ...
    'ContentType','vector');
savefig(simple_fig,strrep(simple_file,'.png','.fig'));
close(simple_fig);

Summary = table(gains,StartErrorMean,TerminalErrorMean,TerminalErrorMax, ...
    TerminalBelowOnePercent,TriggerGapMean,TriggerGapMedian, ...
    TriggerRelativeGapMedian,BroadcastsPerAgent,files, ...
    'VariableNames',{'ConsensusStep','StartErrorMean', ...
    'TerminalErrorMean','TerminalErrorMax','TerminalBelowOnePercent', ...
    'TriggerGapMean','TriggerGapMedian','TriggerRelativeGapMedian', ...
    'BroadcastsPerAgent','ResultFile'});
writetable(Summary,fullfile(result_folder,'gain_tradeoff_summary.csv'));
disp(Summary(:,1:end-1));
fprintf('Saved reference trade-off plot: %s\n',plot_file);
fprintf('Saved selected-window reference plot: %s\n',simple_file);
end

function restore_environment_local(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
