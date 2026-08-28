function Summary = run_ip_ac_periodic_comparison(do_run)
%RUN_IP_AC_PERIODIC_COMPARISON Compare Seyboth and periodic AC triggers.

if nargin < 1 || isempty(do_run), do_run = true; end
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result','ip_ac_periodic_m600');
if ~isfolder(result_folder), mkdir(result_folder); end
point_file = fullfile(repo_root,'result','ip_reachable_m600', ...
    'reachable_uniform_M600.mat');

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_IP_PROJECTION_DIAGNOSTICS','CONTROL_IP_INDUCING_POINT_FILE', ...
    'CONTROL_AC_ITERATION_POLICY','CONTROL_AC_FIXED_ITERATIONS', ...
    'CONTROL_AC_BROADCAST_TRIGGER','CONTROL_AC_PERIODIC_SIGMA'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment_local( ...
    env_names,old_values));

setenv(env_names{1},'30');
setenv(env_names{2},'600');
setenv(env_names{3},'0');
setenv(env_names{4},point_file);
setenv(env_names{5},'fixed');
setenv(env_names{6},'10');
setenv(env_names{8},'0.2');

trigger_modes = {'petc_reset','petc_cache'};
petc_files = strings(2,1);
for case_i = 1:2
    setenv(env_names{7},trigger_modes{case_i});
    name = sprintf('poe_ac_%s_R10_sigma_0p2_T30',trigger_modes{case_i});
    petc_files(case_i) = fullfile(result_folder,[name '.mat']);
    if do_run || ~isfile(petc_files(case_i))
        run_simulation_inducing_point('poe_ac',result_folder,name, ...
            false,[],[],[],42);
    end
end

seyboth_file = fullfile(repo_root,'result','ip_ac_seyboth_m600', ...
    'poe_ac_seyboth_R10_c0_0p5_c1_1_lam_5_T30.mat');
files = [string(seyboth_file); petc_files];
Method = ["Seyboth-AC10"; "PETC-reset-AC10"; "PETC-cache-AC10"];

MeanTrackingError = nan(3,1);
FinalTrackingError = nan(3,1);
MeanPredictionError = nan(3,1);
BroadcastsPerAgent = nan(3,1);
FinalConsensusDisagreement = nan(3,1);
AgentBroadcastCounts = nan(3,6);
time_history = cell(3,1);
tracking_history = cell(3,1);
consensus_history = cell(3,1);
trigger_history = cell(3,1);
update_interval = nan(3,1);

for case_i = 1:3
    d = load(files(case_i),'TrackingError_vector', ...
        'prediction_error_norm_vector','ac_broadcasts_per_agent', ...
        'ac_event_count_per_agent','ac_iteration_count_set', ...
        'ac_consensus_disagreement_set','ac_inner_trigger_instance_set', ...
        't_set','ACOnlineUpdateInterval');
    MeanTrackingError(case_i) = mean(d.TrackingError_vector,'omitnan');
    FinalTrackingError(case_i) = d.TrackingError_vector(end);
    MeanPredictionError(case_i) = ...
        mean(d.prediction_error_norm_vector,'omitnan');
    BroadcastsPerAgent(case_i) = d.ac_broadcasts_per_agent;
    AgentBroadcastCounts(case_i,:) = d.ac_event_count_per_agent(:).';
    last_update = find(d.ac_iteration_count_set > 0,1,'last');
    FinalConsensusDisagreement(case_i) = ...
        d.ac_consensus_disagreement_set(last_update);
    time_history{case_i} = d.t_set;
    tracking_history{case_i} = d.TrackingError_vector;
    consensus_history{case_i} = d.ac_consensus_disagreement_set;
    trigger_history{case_i} = ...
        d.ac_inner_trigger_instance_set(:,1:10,:);
    update_interval(case_i) = d.ACOnlineUpdateInterval;
end

star_fig = figure('Visible','off','Color','w', ...
    'Position',[80 50 1150 850]);
tl = tiledlayout(star_fig,3,1,'TileSpacing','compact','Padding','compact');
colors = lines(6);
for case_i = 1:3
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    trigger_data = trigger_history{case_i};
    t_set = time_history{case_i};
    for agent_i = 1:6
        [round_idx,time_idx] = find(squeeze( ...
            trigger_data(agent_i,:,:)));
        event_time = t_set(time_idx);
        event_time = event_time(:)+(round_idx(:)-1)* ...
            (update_interval(case_i)/10);
        plot(ax,event_time,agent_i*ones(size(event_time)),'*', ...
            'Color',colors(agent_i,:),'MarkerSize',4, ...
            'LineStyle','none','DisplayName',sprintf( ...
            'Agent %d, n=%d',agent_i,numel(event_time)));
    end
    xlim(ax,[0 30]); ylim(ax,[0.5 6.5]); yticks(ax,1:6);
    ylabel(ax,'Agent'); title(ax,char(Method(case_i)));
    legend(ax,'Location','eastoutside');
end
xlabel(tl,'Time (s)');
print(star_fig,fullfile(result_folder, ...
    'ip_ac_trigger_star_seyboth_vs_petc.png'),'-dpng','-r160');
close(star_fig);

metric_fig = figure('Visible','off','Color','w', ...
    'Position',[100 80 1050 720]);
mt = tiledlayout(metric_fig,2,1,'TileSpacing','compact','Padding','compact');
styles = {'-','--',':'};
ax1 = nexttile(mt); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
for case_i = 1:3
    plot(ax1,time_history{case_i},tracking_history{case_i}, ...
        'LineWidth',1.5,'LineStyle',styles{case_i}, ...
        'DisplayName',char(Method(case_i)));
end
ylabel(ax1,'Tracking error ||vartheta(t)||_2'); xlim(ax1,[0 30]);
legend(ax1,'Location','northeast'); title(ax1,'Tracking error');

ax2 = nexttile(mt); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
for case_i = 1:3
    c = consensus_history{case_i};
    valid = isfinite(c);
    semilogy(ax2,time_history{case_i}(valid),max(c(valid),eps), ...
        'LineWidth',1.5,'LineStyle',styles{case_i}, ...
        'DisplayName',char(Method(case_i)));
end
xlabel(ax2,'Time (s)'); ylabel(ax2,'End-of-update disagreement');
xlim(ax2,[0 30]); legend(ax2,'Location','best');
title(ax2,'Consensus disagreement after 10 rounds');
print(metric_fig,fullfile(result_folder, ...
    'ip_ac_metrics_seyboth_vs_petc.png'),'-dpng','-r160');
close(metric_fig);

Summary = table(Method,MeanTrackingError,FinalTrackingError, ...
    MeanPredictionError,BroadcastsPerAgent,FinalConsensusDisagreement);
writetable(Summary,fullfile(result_folder, ...
    'ip_ac_seyboth_vs_petc_summary.csv'));
writematrix(AgentBroadcastCounts,fullfile(result_folder, ...
    'ip_ac_seyboth_vs_petc_agent_counts.csv'));
disp(Summary);
end

function restore_environment_local(names,values)
for i = 1:numel(names), setenv(names{i},values{i}); end
end
