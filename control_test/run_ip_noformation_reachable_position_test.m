function Summary = run_ip_noformation_reachable_position_test(do_run)
%RUN_IP_NOFORMATION_REACHABLE_POSITION_TEST Test known leader-only domain.

if nargin < 1 || isempty(do_run), do_run = true; end
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result', ...
    'ip_noformation_reachable_position_test');
if ~isfolder(result_folder), mkdir(result_folder); end

% Without formation, the prescribed position is the unit leader circle.
% Add a 10% safety margin without using any simulated state trajectory.
position_domain_scale = 1.10;
point_file = fullfile(result_folder, ...
    'inducing_points_position100_velocity6_leader_reachable.mat');
generate_reference_velocity_inducing_points(point_file,100,6, ...
    position_domain_scale);

env_names = {'CONTROL_SIM_END_TIME','CONTROL_ONLINE_AGENT_POLICY', ...
    'CONTROL_IP_NUM_INDUCING_POINTS','CONTROL_IP_INDUCING_POINT_FILE', ...
    'CONTROL_IP_POSITION_LENGTH_SCALE','CONTROL_IP_VELOCITY_LENGTH_SCALE', ...
    'CONTROL_IP_RECONSTRUCTION_SIGMA_N','CONTROL_IP_PROJECTION_DIAGNOSTICS', ...
    'CONTROL_DAC_STEPS_PER_PHYSICAL_STEP','CONTROL_DAC_TRIGGER_EPSILON', ...
    'CONTROL_AGGREGATION_VARIANCE_FLOOR','CONTROL_AGGREGATION_PRECISION_FLOOR', ...
    'CONTROL_GPOE_BETA_MAX','CONTROL_RBCM_BETA_MAX','CONTROL_BCM_PRIOR_SCALE'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment(env_names,old_values)); %#ok<NASGU>
env_values = {'20','all_agents','600',point_file,'0.5','1.5', ...
    '0.000001','0','1','0.005','0.001','0.0001','10','0.25','0.25'};
for env_i = 1:numel(env_names), setenv(env_names{env_i},env_values{env_i}); end

result_file = fullfile(result_folder, ...
    'poe_ip_dac_leader_reachable.mat');
if do_run
    run_simulation_inducing_point('poe',result_folder, ...
        'poe_ip_dac_leader_reachable',false,0.4,0.1,[],1);
end

baseline_file = fullfile(repo_root,'result', ...
    'ip_reference_velocity_validation','poe_ip_dac.mat');
labels = {'Full [-1.5,1.5]';'Leader range + 10% margin'};
files = {baseline_file;result_file};
MaxTrackingErrorAfter10s = nan(2,1);
MeanPredictionError = nan(2,1);
LatePChange = nan(2,1);
BroadcastsPerAgent = nan(2,1);
for case_i = 1:2
    data = load(files{case_i});
    mask = data.t_set >= 10;
    MaxTrackingErrorAfter10s(case_i) = max( ...
        data.TrackingError_vector(mask),[],'omitnan');
    MeanPredictionError(case_i) = mean( ...
        data.prediction_error_norm_vector,'omitnan');
    LatePChange(case_i) = data.LatePChange;
    BroadcastsPerAgent(case_i) = mean( ...
        data.dac_broadcast_event_count_per_agent,'omitnan');
end
SamplingDomain = labels;
Summary = table(SamplingDomain,MaxTrackingErrorAfter10s, ...
    MeanPredictionError,LatePChange,BroadcastsPerAgent);
writetable(Summary,fullfile(result_folder,'summary.csv'));
save(fullfile(result_folder,'summary.mat'),'Summary');
disp(Summary);
end

function restore_environment(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
