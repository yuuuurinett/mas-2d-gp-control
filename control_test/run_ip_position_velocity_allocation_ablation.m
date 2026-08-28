function Summary = run_ip_position_velocity_allocation_ablation(do_run)
%RUN_IP_POSITION_VELOCITY_ALLOCATION_ABLATION Compare three M=600 layouts.
% The unknown dynamics is strongly position dominated, so the comparison
% moves inducing-point budget from velocity pairs to position pairs while
% leaving every simulation/consensus parameter unchanged.

if nargin < 1 || isempty(do_run), do_run = true; end

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result', ...
    'ip_position_velocity_allocation_ablation');
if ~isfolder(result_folder), mkdir(result_folder); end

position_counts = [100 120 150];
velocity_counts = [6 5 4];

env_names = { ...
    'CONTROL_SIM_END_TIME','CONTROL_ONLINE_AGENT_POLICY', ...
    'CONTROL_IP_NUM_INDUCING_POINTS','CONTROL_IP_INDUCING_POINT_FILE', ...
    'CONTROL_IP_POSITION_LENGTH_SCALE','CONTROL_IP_VELOCITY_LENGTH_SCALE', ...
    'CONTROL_IP_RECONSTRUCTION_SIGMA_N','CONTROL_IP_PROJECTION_DIAGNOSTICS', ...
    'CONTROL_DAC_STEPS_PER_PHYSICAL_STEP','CONTROL_DAC_TRIGGER_EPSILON', ...
    'CONTROL_AGGREGATION_VARIANCE_FLOOR','CONTROL_AGGREGATION_PRECISION_FLOOR', ...
    'CONTROL_GPOE_BETA_MAX','CONTROL_RBCM_BETA_MAX', ...
    'CONTROL_BCM_PRIOR_SCALE'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment(env_names,old_values)); %#ok<NASGU>

env_values = {'20','all_agents','600','', '0.5','1.5','0.000001','0', ...
    '1','0.005','0.001','0.0001','10','0.25','0.25'};
for env_i = 1:numel(env_names)
    setenv(env_names{env_i},env_values{env_i});
end

if do_run
    for allocation_i = 1:numel(position_counts)
        n_position = position_counts(allocation_i);
        n_velocity = velocity_counts(allocation_i);
        inducing_file = fullfile(result_folder,sprintf( ...
            'inducing_points_position%d_velocity%d_reference_seed42.mat', ...
            n_position,n_velocity));
        generate_reference_velocity_inducing_points(inducing_file, ...
            n_position,n_velocity);
        setenv('CONTROL_IP_INDUCING_POINT_FILE',inducing_file);

        result_name = sprintf('poe_ip_dac_position%d_velocity%d', ...
            n_position,n_velocity);
        fprintf('\n[Allocation %d/%d] %d position x %d velocity = 600\n', ...
            allocation_i,numel(position_counts),n_position,n_velocity);
        run_simulation_inducing_point('poe',result_folder,result_name, ...
            true,0.4,0.1,[],1);
    end
end

PositionPairs = position_counts(:);
VelocityPairs = velocity_counts(:);
MaxTrackingErrorAfter10s = nan(numel(position_counts),1);
MeanPredictionError = nan(numel(position_counts),1);
MeanDirectPoEError = nan(numel(position_counts),1);
MeanIPCENError = nan(numel(position_counts),1);
BroadcastsPerAgent = nan(numel(position_counts),1);
ElapsedSeconds = nan(numel(position_counts),1);

for allocation_i = 1:numel(position_counts)
    result_file = fullfile(result_folder,sprintf( ...
        'poe_ip_dac_position%d_velocity%d.mat', ...
        position_counts(allocation_i),velocity_counts(allocation_i)));
    data = load(result_file);
    evaluation_mask = data.t_set >= 10;
    MaxTrackingErrorAfter10s(allocation_i) = max( ...
        data.TrackingError_vector(evaluation_mask),[],'omitnan');
    MeanPredictionError(allocation_i) = mean( ...
        data.prediction_error_norm_vector,'omitnan');
    MeanDirectPoEError(allocation_i) = mean( ...
        data.direct_prediction_error_norm_vector,'omitnan');
    MeanIPCENError(allocation_i) = mean( ...
        data.ip_cen_prediction_error_norm_vector,'omitnan');
    BroadcastsPerAgent(allocation_i) = mean( ...
        data.dac_broadcast_event_count_per_agent,'omitnan');
    ElapsedSeconds(allocation_i) = data.elapsed_time;
end

Summary = table(PositionPairs,VelocityPairs,MaxTrackingErrorAfter10s, ...
    MeanPredictionError,MeanDirectPoEError,MeanIPCENError, ...
    BroadcastsPerAgent,ElapsedSeconds);
writetable(Summary,fullfile(result_folder,'summary.csv'));
save(fullfile(result_folder,'summary.mat'),'Summary');
disp(Summary);
end

function restore_environment(names,values)
for idx = 1:numel(names)
    setenv(names{idx},values{idx});
end
end
