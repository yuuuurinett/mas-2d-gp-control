function Summary = run_ip_ac_seyboth_threshold_sweep(do_run)
%RUN_IP_AC_SEYBOTH_THRESHOLD_SWEEP Short AC-10 communication sweep.

if nargin < 1 || isempty(do_run), do_run = true; end
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result','ip_ac_seyboth_threshold_sweep');
if ~isfolder(result_folder), mkdir(result_folder); end
point_file = fullfile(repo_root,'result','ip_reachable_m600', ...
    'reachable_uniform_M600.mat');

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_IP_PROJECTION_DIAGNOSTICS','CONTROL_IP_INDUCING_POINT_FILE', ...
    'CONTROL_AC_ITERATION_POLICY','CONTROL_AC_FIXED_ITERATIONS', ...
    'CONTROL_AC_BROADCAST_TRIGGER','CONTROL_AC_SEYBOTH_C0', ...
    'CONTROL_AC_SEYBOTH_C1','CONTROL_AC_SEYBOTH_LAMBDA', ...
    'CONTROL_AC_TRIGGER_DIAGNOSTICS'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment_local( ...
    env_names,old_values));

setenv(env_names{1},'10');
setenv(env_names{2},'600');
setenv(env_names{3},'0');
setenv(env_names{4},point_file);
setenv(env_names{5},'fixed');
setenv(env_names{6},'10');
setenv(env_names{7},'seyboth');
setenv(env_names{10},'5');
% Save the per-agent, per-inner-round detector measure and threshold so the
% selected T10 case can be plotted in the same form as Fig. 8.
setenv(env_names{11},'1');

c0_values = [0.5;0.75;1.0;1.25];
c1_values = [1.0;1.5;2.0;2.5];
case_count = numel(c0_values);
files = strings(case_count,1);

for case_i = 1:case_count
    setenv(env_names{8},num2str(c0_values(case_i),'%.12g'));
    setenv(env_names{9},num2str(c1_values(case_i),'%.12g'));
    name = sprintf('poe_ac_R10_c0_%s_c1_%s_T10', ...
        number_tag(c0_values(case_i)),number_tag(c1_values(case_i)));
    files(case_i) = fullfile(result_folder,[name '.mat']);
    if do_run || ~isfile(files(case_i))
        run_simulation_inducing_point('poe_ac',result_folder,name, ...
            false,[],[],[],42);
    end
end

MeanTrackingError = nan(case_count,1);
FinalTrackingError = nan(case_count,1);
MeanPredictionError = nan(case_count,1);
MaxPredictionError = nan(case_count,1);
BroadcastsPerAgent = nan(case_count,1);
FinalConsensusDisagreement = nan(case_count,1);
NegativeDenominatorSamples = nan(case_count,1);

for case_i = 1:case_count
    d = load(files(case_i),'TrackingError_vector', ...
        'prediction_error_norm_vector','ac_broadcasts_per_agent', ...
        'ac_iteration_count_set','ac_consensus_disagreement_set', ...
        'aggregation_negative_den_count_vector');
    MeanTrackingError(case_i) = mean(d.TrackingError_vector,'omitnan');
    FinalTrackingError(case_i) = d.TrackingError_vector(end);
    MeanPredictionError(case_i) = ...
        mean(d.prediction_error_norm_vector,'omitnan');
    MaxPredictionError(case_i) = ...
        max(d.prediction_error_norm_vector,[],'omitnan');
    BroadcastsPerAgent(case_i) = d.ac_broadcasts_per_agent;
    last_update = find(d.ac_iteration_count_set > 0,1,'last');
    FinalConsensusDisagreement(case_i) = ...
        d.ac_consensus_disagreement_set(last_update);
    NegativeDenominatorSamples(case_i) = ...
        sum(d.aggregation_negative_den_count_vector,'omitnan');
end

Summary = table(c0_values,c1_values,MeanTrackingError, ...
    FinalTrackingError,MeanPredictionError,MaxPredictionError, ...
    BroadcastsPerAgent,FinalConsensusDisagreement, ...
    NegativeDenominatorSamples, ...
    'VariableNames',{'C0','C1','MeanTrackingError', ...
    'FinalTrackingError','MeanPredictionError','MaxPredictionError', ...
    'BroadcastsPerAgent','FinalConsensusDisagreement', ...
    'NegativeDenominatorSamples'});
writetable(Summary,fullfile(result_folder, ...
    'seyboth_threshold_sweep_T10.csv'));
disp(Summary);
end

function tag = number_tag(value)
tag = strrep(num2str(value,'%.3g'),'.','p');
end

function restore_environment_local(names,values)
for i = 1:numel(names), setenv(names{i},values{i}); end
end
