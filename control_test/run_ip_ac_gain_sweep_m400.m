function Summary = run_ip_ac_gain_sweep_m400(do_run)
%RUN_IP_AC_GAIN_SWEEP_M400 Fixed-R10 AC gain comparison.
% All cases use M=400, sigma=0.3, 300 online samples per agent, and T=30 s.

if nargin < 1 || isempty(do_run), do_run = true; end
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result','ip_ac_gain_sweep_m400_T30');
if ~isfolder(result_folder), mkdir(result_folder); end

gains = [0.1;0.2;0.3;0.4];
case_count = numel(gains);
files = strings(case_count,1);

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_IP_INDUCING_POINT_FILE','CONTROL_AC_ITERATION_POLICY', ...
    'CONTROL_AC_FIXED_ITERATIONS','CONTROL_AC_BROADCAST_TRIGGER', ...
    'CONTROL_AC_PERIODIC_SIGMA','CONTROL_AC_TRIGGER_DIAGNOSTICS', ...
    'CONTROL_ONLINE_AGENT_POLICY','CONTROL_AC_CONSENSUS_STEP'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment_local(env_names,old_values)); %#ok<NASGU>

setenv(env_names{1},'30');
setenv(env_names{2},'400');
setenv(env_names{3},'');
setenv(env_names{4},'fixed');
setenv(env_names{5},'10');
setenv(env_names{6},'petc_cache');
setenv(env_names{7},'0.3');
setenv(env_names{8},'0');
setenv(env_names{9},'all_agents');

for case_i = 1:case_count
    setenv(env_names{10},num2str(gains(case_i),'%.12g'));
    name = sprintf('poe_ac_gain_%s_R10_M400_online300_per_agent_T30', ...
        number_tag(gains(case_i)));
    files(case_i) = fullfile(result_folder,[name '.mat']);
    if do_run || ~isfile(files(case_i))
        run_simulation_inducing_point('poe_ac',result_folder,name, ...
            false,[],[],[],42);
    end
end

Round10ErrorMean = nan(case_count,1);
Round10ErrorMax = nan(case_count,1);
Round10ErrorAfter10sMean = nan(case_count,1);
Round10ErrorAfter10sMax = nan(case_count,1);
MaxTrackingErrorAfter10s = nan(case_count,1);
MeanPredictionError = nan(case_count,1);
BroadcastsPerAgent = nan(case_count,1);
TotalBroadcasts = nan(case_count,1);

for case_i = 1:case_count
    d = load(files(case_i),'t_set','TrackingError_vector', ...
        'prediction_error_norm_vector','ac_event_count_per_agent', ...
        'ac_iteration_count_set','ac_tracking_error_history_set');
    update_idx = find(d.ac_iteration_count_set > 0);
    round10_error = nan(size(update_idx));
    for idx = 1:numel(update_idx)
        history = d.ac_tracking_error_history_set{update_idx(idx)};
        if ~isempty(history), round10_error(idx) = history(end); end
    end
    after10_update = d.t_set(update_idx) >= 10;
    Round10ErrorMean(case_i) = mean(round10_error,'omitnan');
    Round10ErrorMax(case_i) = max(round10_error,[],'omitnan');
    Round10ErrorAfter10sMean(case_i) = ...
        mean(round10_error(after10_update),'omitnan');
    Round10ErrorAfter10sMax(case_i) = ...
        max(round10_error(after10_update),[],'omitnan');
    physical_after10 = d.t_set >= 10;
    MaxTrackingErrorAfter10s(case_i) = max( ...
        d.TrackingError_vector(physical_after10),[],'omitnan');
    MeanPredictionError(case_i) = ...
        mean(d.prediction_error_norm_vector,'omitnan');
    TotalBroadcasts(case_i) = sum(d.ac_event_count_per_agent);
    BroadcastsPerAgent(case_i) = mean(d.ac_event_count_per_agent);
end

% For the current ring graph, lambda_2=1 and lambda_max=4.
StaticContractionFactor = max(abs([1-gains,1-4*gains]),[],2);
Summary = table(gains,StaticContractionFactor,Round10ErrorMean, ...
    Round10ErrorMax,Round10ErrorAfter10sMean,Round10ErrorAfter10sMax, ...
    MaxTrackingErrorAfter10s,MeanPredictionError,BroadcastsPerAgent, ...
    TotalBroadcasts,files, ...
    'VariableNames',{'ConsensusGain','StaticContractionFactor', ...
    'Round10ErrorMean','Round10ErrorMax','Round10ErrorAfter10sMean', ...
    'Round10ErrorAfter10sMax','MaxTrackingErrorAfter10s', ...
    'MeanPredictionError','BroadcastsPerAgent','TotalBroadcasts', ...
    'ResultFile'});
writetable(Summary,fullfile(result_folder,'ac_gain_sweep_summary.csv'));
disp(Summary(:,1:end-1));
end

function tag = number_tag(value)
tag = strrep(num2str(value,'%.3g'),'.','p');
end

function restore_environment_local(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
