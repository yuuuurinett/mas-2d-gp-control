function Summary = run_ip_ac_meng_sigma_sweep(do_run,sigma_values)
%RUN_IP_AC_MENG_SIGMA_SWEEP AC-10 communication sweep over Meng et al.
%(2015) periodic cached-broadcast sigma at T=30s.

if nargin < 1 || isempty(do_run), do_run = true; end
if nargin < 2 || isempty(sigma_values)
    sigma_values = [0.20;0.25;0.30;0.35];
else
    sigma_values = sigma_values(:);
end
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result', ...
    'ip_ac_meng_sigma_sweep_m400_online300_per_agent');
if ~isfolder(result_folder), mkdir(result_folder); end

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_IP_PROJECTION_DIAGNOSTICS','CONTROL_IP_INDUCING_POINT_FILE', ...
    'CONTROL_AC_FIXED_ITERATIONS','CONTROL_AC_PERIODIC_SIGMA', ...
    'CONTROL_AC_TRIGGER_DIAGNOSTICS','CONTROL_ONLINE_AGENT_POLICY'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment_local( ...
    env_names,old_values));

setenv(env_names{1},'30');
setenv(env_names{2},'400');
setenv(env_names{3},'0');
setenv(env_names{4},'');                % uniform random inducing points
setenv(env_names{5},'10');
setenv(env_names{7},'1');               % save LHS/RHS for the detector plot
setenv(env_names{8},'all_agents');       % 300 online samples per agent

case_count = numel(sigma_values);
files = strings(case_count,1);

for case_i = 1:case_count
    setenv(env_names{6},num2str(sigma_values(case_i),'%.12g'));
    name = sprintf('poe_ac_petc_R10_M400_online300_per_agent_sigma_%s_T30', ...
        number_tag(sigma_values(case_i)));
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

Summary = table(sigma_values,MeanTrackingError, ...
    FinalTrackingError,MeanPredictionError,MaxPredictionError, ...
    BroadcastsPerAgent,FinalConsensusDisagreement, ...
    NegativeDenominatorSamples, ...
    'VariableNames',{'Sigma','MeanTrackingError', ...
    'FinalTrackingError','MeanPredictionError','MaxPredictionError', ...
    'BroadcastsPerAgent','FinalConsensusDisagreement', ...
    'NegativeDenominatorSamples'});
writetable(Summary,fullfile(result_folder, ...
    'meng_sigma_sweep_T30.csv'));
disp(Summary);

% Generate the advisor-facing detector figure for sigma=0.3 from the same
% result file and the exact LHS/RHS values used by the simulation.
sigma03_idx = find(abs(sigma_values-0.3) < 1e-12,1,'first');
if ~isempty(sigma03_idx)
    for agent_i = 1:6
        plot_ip_ac_advisor_detector(files(sigma03_idx),agent_i,10.0);
    end
end
end

function tag = number_tag(value)
tag = strrep(num2str(value,'%.3g'),'.','p');
end

function restore_environment_local(names,values)
for i = 1:numel(names), setenv(names{i},values{i}); end
end
