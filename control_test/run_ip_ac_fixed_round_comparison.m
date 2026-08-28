function Summary = run_ip_ac_fixed_round_comparison(simulation_end_time,reuse_existing)
%RUN_IP_AC_FIXED_ROUND_COMPARISON Compare restarted IP-AC at 10/20 rounds.

if nargin < 1 || isempty(simulation_end_time), simulation_end_time = 30; end
if nargin < 2 || isempty(reuse_existing), reuse_existing = false; end
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result','ip_reachable_m600');
point_file = fullfile(result_folder,'reachable_uniform_M600.mat');
if ~isfile(point_file)
    error('Inducing-point file not found: %s',point_file);
end

round_values = [10 20];
result_files = strings(numel(round_values),1);
env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_IP_PROJECTION_DIAGNOSTICS','CONTROL_IP_INDUCING_POINT_FILE', ...
    'CONTROL_DAC_KAPPA_P','CONTROL_AC_TRIGGER_SIGMA', ...
    'CONTROL_AC_CONSENSUS_STEP','CONTROL_AC_ITERATION_POLICY', ...
    'CONTROL_AC_FIXED_ITERATIONS'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment(env_names,old_values)); %#ok<NASGU>

setenv(env_names{1},num2str(simulation_end_time,'%.15g'));
setenv(env_names{2},'600');
setenv(env_names{3},'0');
setenv(env_names{4},point_file);
setenv(env_names{5},'1');
setenv(env_names{6},'0.9');
setenv(env_names{7},'0.1');
setenv(env_names{8},'fixed');

for case_idx = 1:numel(round_values)
    rounds = round_values(case_idx);
    setenv(env_names{9},num2str(rounds));
    base_name = sprintf('poe_ac_M600_reachable_fixed_R%d_T%g', ...
        rounds,simulation_end_time);
    result_files(case_idx) = fullfile(result_folder,[base_name '.mat']);
    if ~(reuse_existing && isfile(result_files(case_idx)))
        run_simulation_inducing_point('poe_ac',result_folder,base_name, ...
            false,[],[],[],42);
    end
end

MeanFinalConsensusError = nan(2,1);
MedianFinalConsensusError = nan(2,1);
P90FinalConsensusError = nan(2,1);
MaxFinalConsensusError = nan(2,1);
MeanBroadcastsPerAgent = nan(2,1);
MeanTrackingError = nan(2,1);
MeanPredictionError = nan(2,1);

fig = figure('Color','w','Position',[100 100 850 520]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on');
ax.Toolbar.Visible = 'off';
colors = lines(2);
for case_idx = 1:2
    data = load(result_files(case_idx));
    histories = data.ac_tracking_error_history_set;
    histories = histories(~cellfun(@isempty,histories));
    history_matrix = vertcat(histories{:});
    final_errors = history_matrix(:,end);
    MeanFinalConsensusError(case_idx) = mean(final_errors,'omitnan');
    MedianFinalConsensusError(case_idx) = median(final_errors,'omitnan');
    sorted_errors = sort(final_errors(isfinite(final_errors)));
    p90_idx = max(1,ceil(0.9*numel(sorted_errors)));
    P90FinalConsensusError(case_idx) = sorted_errors(p90_idx);
    MaxFinalConsensusError(case_idx) = max(final_errors,[],'omitnan');
    MeanBroadcastsPerAgent(case_idx) = data.ac_broadcasts_per_agent;
    MeanTrackingError(case_idx) = mean(data.TrackingError_vector,'omitnan');
    MeanPredictionError(case_idx) = mean(data.prediction_error_norm_vector,'omitnan');
    plot(ax,0:round_values(case_idx),mean(history_matrix,1,'omitnan'), ...
        'LineWidth',1.8,'Color',colors(case_idx,:), ...
        'DisplayName',sprintf('AC-%d',round_values(case_idx)));
end
xlabel(ax,'Consensus round');
ylabel(ax,'Mean RMS error to exact P average');
title(ax,'Restarted IP-AC convergence after each P update');
legend(ax,'Location','northeast');
set(ax,'YScale','log');
exportgraphics(fig,fullfile(result_folder, ...
    sprintf('ip_ac_fixed_R10_R20_T%g_convergence.png',simulation_end_time)), ...
    'Resolution',180);

Rounds = round_values(:);
Summary = table(Rounds,MeanFinalConsensusError,MedianFinalConsensusError, ...
    P90FinalConsensusError,MaxFinalConsensusError,MeanBroadcastsPerAgent, ...
    MeanTrackingError,MeanPredictionError);
writetable(Summary,fullfile(result_folder, ...
    sprintf('ip_ac_fixed_R10_R20_T%g_summary.csv',simulation_end_time)));
disp(Summary);
end

function restore_environment(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
