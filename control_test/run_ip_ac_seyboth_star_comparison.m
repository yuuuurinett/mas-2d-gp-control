function Summary = run_ip_ac_seyboth_star_comparison(do_run)
%RUN_IP_AC_SEYBOTH_STAR_COMPARISON AC-10/20 with local Seyboth broadcasts.

if nargin < 1 || isempty(do_run), do_run = true; end
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result','ip_ac_seyboth_m600');
if ~isfolder(result_folder), mkdir(result_folder); end
point_file = fullfile(repo_root,'result','ip_reachable_m600', ...
    'reachable_uniform_M600.mat');

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_IP_PROJECTION_DIAGNOSTICS','CONTROL_IP_INDUCING_POINT_FILE', ...
    'CONTROL_AC_ITERATION_POLICY','CONTROL_AC_FIXED_ITERATIONS', ...
    'CONTROL_AC_SEYBOTH_C0','CONTROL_AC_SEYBOTH_C1', ...
    'CONTROL_AC_SEYBOTH_LAMBDA'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment_local( ...
    env_names,old_values));

setenv(env_names{1},'30');
setenv(env_names{2},'600');
setenv(env_names{3},'0');
setenv(env_names{4},point_file);
setenv(env_names{5},'fixed');
setenv(env_names{7},'0.5');
setenv(env_names{8},'1');
setenv(env_names{9},'5');

round_values = [10 20];
files = strings(2,1);
for case_i = 1:2
    rounds = round_values(case_i);
    setenv(env_names{6},num2str(rounds));
    name = sprintf('poe_ac_seyboth_R%d_c0_0p5_c1_1_lam_5_T30',rounds);
    files(case_i) = fullfile(result_folder,[name '.mat']);
    if do_run || ~isfile(files(case_i))
        run_simulation_inducing_point('poe_ac',result_folder,name, ...
            false,[],[],[],42);
    end
end

Method = strings(2,1);
MeanTrackingError = nan(2,1);
FinalTrackingError = nan(2,1);
MeanPredictionError = nan(2,1);
BroadcastsPerAgent = nan(2,1);
FinalConsensusDisagreement = nan(2,1);
time_history = cell(2,1);
tracking_history = cell(2,1);

star_fig = figure('Visible','off','Color','w','Position',[80 80 1120 700]);
tl = tiledlayout(star_fig,2,1,'TileSpacing','compact','Padding','compact');
colors = lines(6);
for case_i = 1:2
    d = load(files(case_i),'TrackingError_vector', ...
        'prediction_error_norm_vector','ac_broadcasts_per_agent', ...
        'ac_iteration_count_set','ac_consensus_disagreement_set', ...
        'ac_inner_trigger_instance_set','t_set', ...
        'ACOnlineUpdateInterval');
    rounds = round_values(case_i);
    Method(case_i) = sprintf('IP-AC-%d Seyboth',rounds);
    MeanTrackingError(case_i) = mean(d.TrackingError_vector,'omitnan');
    FinalTrackingError(case_i) = d.TrackingError_vector(end);
    MeanPredictionError(case_i) = ...
        mean(d.prediction_error_norm_vector,'omitnan');
    BroadcastsPerAgent(case_i) = d.ac_broadcasts_per_agent;
    last_update = find(d.ac_iteration_count_set > 0,1,'last');
    FinalConsensusDisagreement(case_i) = ...
        d.ac_consensus_disagreement_set(last_update);
    time_history{case_i} = d.t_set;
    tracking_history{case_i} = d.TrackingError_vector;

    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    ax.Toolbar.Visible = 'off';
    trigger_data = d.ac_inner_trigger_instance_set(:,1:rounds,:);
    d = rmfield(d,'ac_inner_trigger_instance_set');
    counts = zeros(6,1);
    for agent_i = 1:6
        [round_idx,time_idx] = find(squeeze(trigger_data(agent_i,:,:)));
        event_time = d.t_set(time_idx);
        event_time = event_time(:) + (round_idx(:)-1)* ...
            (d.ACOnlineUpdateInterval/rounds);
        counts(agent_i) = numel(event_time);
        plot(ax,event_time,agent_i*ones(size(event_time)),'*', ...
            'Color',colors(agent_i,:),'MarkerSize',4, ...
            'LineStyle','none','DisplayName',sprintf( ...
            'Agent %d, n=%d',agent_i,counts(agent_i)));
    end
    xlim(ax,[0 30]); ylim(ax,[0.5 6.5]); yticks(ax,1:6);
    ylabel(ax,'Agent');
    title(ax,sprintf(['IP-AC-%d Seyboth broadcast events: ', ...
        'c_0=0.5, c_1=1, lambda=5'],rounds));
    legend(ax,'Location','eastoutside');
end
xlabel(tl,'Time (s)');
print(star_fig,fullfile(result_folder, ...
    'ip_ac_seyboth_star_triggers_R10_R20.png'),'-dpng','-r160');
close(star_fig);

track_fig = figure('Visible','off','Color','w','Position',[100 100 900 520]);
ax = axes(track_fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
styles = {'-','--'};
for case_i = 1:2
    semilogy(ax,time_history{case_i}, ...
        max(tracking_history{case_i},eps), ...
        'LineWidth',1.6,'LineStyle',styles{case_i}, ...
        'DisplayName',char(Method(case_i)));
end
xlabel(ax,'Time (s)'); ylabel(ax,'Overall tracking error ||vartheta(t)||_2');
title(ax,'IP-AC with local Seyboth broadcast trigger');
xlim(ax,[0 30]); legend(ax,'Location','northeast');
print(track_fig,fullfile(result_folder, ...
    'ip_ac_seyboth_tracking_R10_R20.png'),'-dpng','-r160');
close(track_fig);

Summary = table(Method,MeanTrackingError,FinalTrackingError, ...
    MeanPredictionError,BroadcastsPerAgent,FinalConsensusDisagreement);
writetable(Summary,fullfile(result_folder,'ip_ac_seyboth_summary.csv'));
disp(Summary);
end

function restore_environment_local(names,values)
for i = 1:numel(names), setenv(names{i},values{i}); end
end
