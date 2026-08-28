function Summary = run_ip_local_unknown_disturbance_tradeoff(do_run)
%RUN_IP_LOCAL_UNKNOWN_DISTURBANCE_TRADEOFF Focused 20-s sensitivity test.
% Every case uses identical data-count conventions, inducing points,
% consensus settings and initial seed. Only the plant/label scales change.

if nargin < 1 || isempty(do_run), do_run = true; end
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result', ...
    'ip_local_unknown_disturbance_tradeoff_T20');
if ~isfolder(result_folder), mkdir(result_folder); end

unknown_values = [0.2;0.3;0.4;0.3;0.4;0.6;0.8;1.0];
disturbance_values = [0.1;0.1;0.1;0.5;0.5;0.1;0.1;0.1];
case_count = numel(unknown_values);

point_file = fullfile(repo_root,'result', ...
    'inducing_points_position_grid100_velocity6_full_domain.mat');
env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_IP_INDUCING_POINT_FILE','CONTROL_IP_POSITION_LENGTH_SCALE', ...
    'CONTROL_IP_VELOCITY_LENGTH_SCALE', ...
    'CONTROL_DAC_STEPS_PER_PHYSICAL_STEP', ...
    'CONTROL_DAC_TRIGGER_EPSILON'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment_local( ...
    env_names,old_values)); %#ok<NASGU>
setenv('CONTROL_SIM_END_TIME','20');
setenv('CONTROL_IP_NUM_INDUCING_POINTS','600');
setenv('CONTROL_IP_INDUCING_POINT_FILE',point_file);
setenv('CONTROL_IP_POSITION_LENGTH_SCALE','0.5');
setenv('CONTROL_IP_VELOCITY_LENGTH_SCALE','1.5');
setenv('CONTROL_DAC_STEPS_PER_PHYSICAL_STEP','1');
setenv('CONTROL_DAC_TRIGGER_EPSILON','0.005');

ip_files = strings(case_count,1);
local_files = strings(case_count,1);
for case_i = 1:case_count
    u_tag = strrep(sprintf('%.1f',unknown_values(case_i)),'.','p');
    d_tag = strrep(sprintf('%.1f',disturbance_values(case_i)),'.','p');
    tag = sprintf('u%s_d%s',u_tag,d_tag);
    ip_files(case_i) = fullfile(result_folder, ...
        ['ip_dac_' tag '_seed42.mat']);
    local_files(case_i) = fullfile(result_folder, ...
        ['offline_local_' tag '_seed42.mat']);
    if do_run && ~isfile(ip_files(case_i))
        run_simulation_inducing_point('poe',result_folder, ...
            ['ip_dac_' tag '_seed42'],false, ...
            unknown_values(case_i),disturbance_values(case_i),[],42);
    end
    if do_run && ~isfile(local_files(case_i))
        run_simulation_test_point('local_offline',result_folder, ...
            ['offline_local_' tag '_seed42'],false,42,350, ...
            unknown_values(case_i),disturbance_values(case_i));
    end
end

IPMaxTrackingAfter10 = nan(case_count,1);
LocalMaxTrackingAfter10 = nan(case_count,1);
TrackingImprovementPercent = nan(case_count,1);
IPMeanPredictionAfter10 = nan(case_count,1);
LocalMeanPredictionAfter10 = nan(case_count,1);
IPBroadcastsPerAgent = nan(case_count,1);
IPMinAggregationDenominator = nan(case_count,1);
IPMaxAbsAggregationTarget = nan(case_count,1);

for case_i = 1:case_count
    ip = load(ip_files(case_i));
    local = load(local_files(case_i));
    ip_tracking = sqrt(sum(ip.vartheta_all_set.^2,1));
    local_tracking = sqrt(sum(local.vartheta_all_set.^2,1));
    ip_mask = ip.t_set >= 10;
    local_mask = local.t_set >= 10;
    IPMaxTrackingAfter10(case_i) = max(ip_tracking(ip_mask),[],'omitnan');
    LocalMaxTrackingAfter10(case_i) = ...
        max(local_tracking(local_mask),[],'omitnan');
    TrackingImprovementPercent(case_i) = 100*( ...
        LocalMaxTrackingAfter10(case_i)-IPMaxTrackingAfter10(case_i)) ...
        /LocalMaxTrackingAfter10(case_i);
    IPMeanPredictionAfter10(case_i) = mean( ...
        ip.prediction_error_norm_vector(ip_mask),'omitnan');
    LocalMeanPredictionAfter10(case_i) = mean( ...
        local.prediction_error_norm_vector(local_mask),'omitnan');
    IPBroadcastsPerAgent(case_i) = ip.dac_broadcasts_per_agent;
    IPMinAggregationDenominator(case_i) = min( ...
        ip.aggregation_min_abs_den_vector,[],'omitnan');
    IPMaxAbsAggregationTarget(case_i) = max( ...
        ip.aggregation_max_abs_phi_vector,[],'omitnan');
end

UnknownScale = unknown_values;
DisturbanceScale = disturbance_values;
Summary = table(UnknownScale,DisturbanceScale, ...
    IPMaxTrackingAfter10,LocalMaxTrackingAfter10, ...
    TrackingImprovementPercent,IPMeanPredictionAfter10, ...
    LocalMeanPredictionAfter10,IPBroadcastsPerAgent, ...
    IPMinAggregationDenominator,IPMaxAbsAggregationTarget, ...
    ip_files,local_files);
writetable(Summary,fullfile(result_folder,'tradeoff_summary.csv'));
save(fullfile(result_folder,'tradeoff_summary.mat'),'Summary');
disp(Summary(:,1:10));
end

function restore_environment_local(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
