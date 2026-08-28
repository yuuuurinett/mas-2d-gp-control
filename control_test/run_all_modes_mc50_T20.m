function [MeanTable,StdTable,LongTable] = ...
    run_all_modes_mc50_T20(do_run,seed_subset,do_summarize,use_formation)
%RUN_ALL_MODES_MC50_T20 Final all-mode Monte Carlo comparison.
% Common configuration:
%   T=20 s, metric=max_{t>=10 s}||vartheta(t)||_2
%   online data=one sample/agent/0.1 s, offline=300 samples/agent
%   unknown scale=0.4, disturbance scale=0.1
%   IP=M600 structured Z (100 position pairs x 6 reference-velocity pairs)
%   TP-DAC distributes ten updates over the ten physical 0.01 s steps.
%   IP-AC and TP-AC restart on each 0.1 s reference and use R=10.

if nargin < 1 || isempty(do_run), do_run = true; end
if nargin < 2 || isempty(seed_subset), seed_subset = 1:50; end
if nargin < 3 || isempty(do_summarize), do_summarize = true; end
if nargin < 4 || isempty(use_formation), use_formation = false; end

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_name = 'mc50_all_modes_T20_U0p4_D0p1_M600_structured_R10';
if use_formation
    result_name = [result_name '_formation'];
end
result_folder = fullfile(repo_root,'result',result_name);
if ~isfolder(result_folder), mkdir(result_folder); end

cfg = final_configuration(repo_root,use_formation);
configure_environment(cfg);
methods = {'poe','gpoe','moe','bcm','rbcm'};
mode_names = {'IP_DAC','IP_AC','TP_DAC','TP_AC','CEN','NBR', ...
    'Offline','Exact'};
seed_values = seed_subset(:).';

if do_run
    for seed = seed_subset(:).'
        if ~ismember(seed,seed_values)
            error('Seed must be in 1:50. Received %d.',seed);
        end
        run_one_seed(seed,methods,result_folder,cfg);
    end
end

if ~do_summarize
    MeanTable = table(); StdTable = table(); LongTable = table();
    return;
end

metric = nan(numel(methods),numel(mode_names),numel(seed_values));
communication = nan(numel(methods),numel(mode_names),numel(seed_values));
for seed_idx = 1:numel(seed_values)
    seed = seed_values(seed_idx);
    seed_folder = fullfile(result_folder,sprintf('seed_%03d',seed));
    offline_value = tracking_metric(fullfile(seed_folder,'offline.mat'));
    exact_value = tracking_metric(fullfile(seed_folder,'exact.mat'));
    for method_idx = 1:numel(methods)
        method = methods{method_idx};
        metric(method_idx,1,seed_idx) = tracking_metric( ...
            fullfile(seed_folder,[method '_ip_dac.mat']));
        metric(method_idx,2,seed_idx) = tracking_metric( ...
            fullfile(seed_folder,[method '_ip_ac.mat']));
        metric(method_idx,3,seed_idx) = tracking_metric( ...
            fullfile(seed_folder,[method '_tp_dac.mat']));
        metric(method_idx,4,seed_idx) = tracking_metric( ...
            fullfile(seed_folder,[method '_tp_ac.mat']));
        metric(method_idx,5,seed_idx) = tracking_metric( ...
            fullfile(seed_folder,[method '_cen.mat']));
        metric(method_idx,6,seed_idx) = tracking_metric( ...
            fullfile(seed_folder,[method '_nbr.mat']));
        metric(method_idx,7,seed_idx) = offline_value;
        metric(method_idx,8,seed_idx) = exact_value;
        communication(method_idx,1,seed_idx) = communication_metric( ...
            fullfile(seed_folder,[method '_ip_dac.mat']),'ip_dac');
        communication(method_idx,2,seed_idx) = communication_metric( ...
            fullfile(seed_folder,[method '_ip_ac.mat']),'ip_ac');
        communication(method_idx,3,seed_idx) = communication_metric( ...
            fullfile(seed_folder,[method '_tp_dac.mat']),'tp');
        communication(method_idx,4,seed_idx) = communication_metric( ...
            fullfile(seed_folder,[method '_tp_ac.mat']),'tp');
        communication(method_idx,6,seed_idx) = communication_metric( ...
            fullfile(seed_folder,[method '_nbr.mat']),'nbr');
    end
end

row_names = cellstr(upper(string(methods)));
MeanTable = array2table(mean(metric,3),'VariableNames',mode_names, ...
    'RowNames',row_names);
StdTable = array2table(std(metric,0,3),'VariableNames',mode_names, ...
    'RowNames',row_names);
formatted = strings(size(metric,1),size(metric,2));
for i = 1:size(formatted,1)
    for j = 1:size(formatted,2)
        formatted(i,j) = sprintf('%.5f +/- %.5f', ...
            MeanTable{i,j},StdTable{i,j});
    end
end
MeanStdTable = array2table(formatted,'VariableNames',mode_names, ...
    'RowNames',row_names);

CommunicationMeanTable = array2table(mean(communication,3,'omitnan'), ...
    'VariableNames',mode_names,'RowNames',row_names);
CommunicationStdTable = array2table(std(communication,0,3,'omitnan'), ...
    'VariableNames',mode_names,'RowNames',row_names);
communication_formatted = strings(size(communication,1), ...
    size(communication,2));
for i = 1:size(communication_formatted,1)
    for j = 1:size(communication_formatted,2)
        if isfinite(CommunicationMeanTable{i,j})
            communication_formatted(i,j) = sprintf('%.2f +/- %.2f', ...
                CommunicationMeanTable{i,j},CommunicationStdTable{i,j});
        else
            communication_formatted(i,j) = "--";
        end
    end
end
CommunicationMeanStdTable = array2table(communication_formatted, ...
    'VariableNames',mode_names,'RowNames',row_names);

row_count = numel(methods)*numel(mode_names)*numel(seed_values);
Aggregation = strings(row_count,1); Mode = strings(row_count,1);
Seed = zeros(row_count,1); MaxTrackingErrorAfter10s = nan(row_count,1);
row = 0;
for i = 1:numel(methods)
    for j = 1:numel(mode_names)
        for k = 1:numel(seed_values)
            row = row+1;
            Aggregation(row) = upper(methods{i});
            Mode(row) = mode_names{j}; Seed(row) = seed_values(k);
            MaxTrackingErrorAfter10s(row) = metric(i,j,k);
        end
    end
end
LongTable = table(Aggregation,Mode,Seed,MaxTrackingErrorAfter10s);

writetable(MeanTable,fullfile(result_folder,'tracking_mean.csv'), ...
    'WriteRowNames',true);
writetable(StdTable,fullfile(result_folder,'tracking_std.csv'), ...
    'WriteRowNames',true);
writetable(MeanStdTable,fullfile(result_folder,'tracking_mean_std.csv'), ...
    'WriteRowNames',true);
writetable(LongTable,fullfile(result_folder,'tracking_all_seeds.csv'));
writetable(CommunicationMeanTable,fullfile(result_folder, ...
    'communication_mean.csv'),'WriteRowNames',true);
writetable(CommunicationStdTable,fullfile(result_folder, ...
    'communication_std.csv'),'WriteRowNames',true);
writetable(CommunicationMeanStdTable,fullfile(result_folder, ...
    'communication_mean_std.csv'),'WriteRowNames',true);
write_complete_markdown(result_folder,MeanTable,StdTable, ...
    CommunicationMeanTable,CommunicationStdTable,methods,mode_names, ...
    seed_values,cfg);
save(fullfile(result_folder,'mc50_summary.mat'),'MeanTable','StdTable', ...
    'MeanStdTable','CommunicationMeanTable','CommunicationStdTable', ...
    'CommunicationMeanStdTable','LongTable','metric','communication', ...
    'methods','mode_names','seed_values','cfg');
disp(MeanStdTable);
disp(CommunicationMeanStdTable);
end

function run_one_seed(seed,methods,result_folder,cfg)
seed_folder = fullfile(result_folder,sprintf('seed_%03d',seed));
if ~isfolder(seed_folder), mkdir(seed_folder); end

for method_idx = 1:numel(methods)
    method = methods{method_idx};
    run_if_invalid(fullfile(seed_folder,[method '_ip_dac.mat']), ...
        'ip',cfg,@() run_simulation_inducing_point(method,seed_folder, ...
        [method '_ip_dac'],cfg.use_formation,cfg.unknown_scale, ...
        cfg.disturbance_scale,[],seed));
    run_if_invalid(fullfile(seed_folder,[method '_ip_ac.mat']), ...
        'ip',cfg,@() run_simulation_inducing_point([method '_ac'], ...
        seed_folder,[method '_ip_ac'],cfg.use_formation,cfg.unknown_scale, ...
        cfg.disturbance_scale,[],seed));
    run_if_invalid(fullfile(seed_folder,[method '_tp_dac.mat']), ...
        'tp',cfg,@() run_simulation_test_point(method,seed_folder, ...
        [method '_tp_dac'],cfg.use_formation,seed,0,cfg.unknown_scale, ...
        cfg.disturbance_scale));
    run_if_invalid(fullfile(seed_folder,[method '_tp_ac.mat']), ...
        'tp',cfg,@() run_simulation_test_point([method '_ac'], ...
        seed_folder,[method '_tp_ac'],cfg.use_formation,seed,0, ...
        cfg.unknown_scale,cfg.disturbance_scale));
    run_if_invalid(fullfile(seed_folder,[method '_cen.mat']), ...
        'aggregation',cfg,@() run_simulation_cen(method,seed_folder, ...
        [method '_cen'],cfg.use_formation,seed,0,cfg.unknown_scale, ...
        cfg.disturbance_scale));
    run_if_invalid(fullfile(seed_folder,[method '_nbr.mat']), ...
        'aggregation',cfg,@() run_simulation_nbr(method,seed_folder, ...
        [method '_nbr'],cfg.use_formation,seed,0,cfg.unknown_scale, ...
        cfg.disturbance_scale));
end

run_if_invalid(fullfile(seed_folder,'offline.mat'),'offline',cfg,@() ...
    run_simulation_test_point('local_offline',seed_folder,'offline', ...
    cfg.use_formation,seed,cfg.offline_per_agent,cfg.unknown_scale, ...
    cfg.disturbance_scale));
run_if_invalid(fullfile(seed_folder,'exact.mat'),'common',cfg,@() ...
    run_simulation_test_point('exact',seed_folder,'exact', ...
    cfg.use_formation,seed,0, ...
    cfg.unknown_scale,cfg.disturbance_scale));
fprintf('Final all-mode MC50 seed %d complete.\n',seed);
end

function run_if_invalid(target,kind,cfg,runner)
if isfile(target) && result_is_valid(target,kind,cfg), return; end
evalc('runner();');
if ~result_is_valid(target,kind,cfg)
    error('New result failed configuration validation: %s',target);
end
end

function valid = result_is_valid(target,kind,cfg)
valid = false;
try
    d = load(target);
    valid = isfield(d,'t_set') && abs(d.t_set(end)-cfg.t_end)<1e-9 && ...
        isfield(d,'UnknownScale') && ...
        abs(d.UnknownScale-cfg.unknown_scale)<1e-12 && ...
        isfield(d,'DisturbanceScale') && ...
        abs(d.DisturbanceScale-cfg.disturbance_scale)<1e-12 && ...
        isfield(d,'TrackingError_vector') && ...
        all(isfinite(d.TrackingError_vector)) && ...
        isfield(d,'use_formation') && ...
        logical(d.use_formation)==logical(cfg.use_formation);
    if ~valid, return; end
    switch kind
        case 'ip'
            expected_points = load(cfg.inducing_file, ...
                'InducingPoints_Coordinates');
            valid = d.NumInducingPoints==600 && ...
                contains(d.InducingPointSource, ...
                'inducing_points_position100_velocity6_reference_seed42.mat') && ...
                isfield(d,'InducingPoints_Coordinates') && ...
                isequal(size(d.InducingPoints_Coordinates), ...
                size(expected_points.InducingPoints_Coordinates)) && ...
                max(abs(d.InducingPoints_Coordinates(:)- ...
                expected_points.InducingPoints_Coordinates(:)))<1e-12 && ...
                d.ACFixedIterations==10 && ...
                abs(d.ACConsensusStep-0.2)<1e-12 && ...
                abs(d.ACPeriodicSigma-0.3)<1e-12 && ...
                isfield(d,'ACForceFullBroadcast') && ...
                ~d.ACForceFullBroadcast;
            if valid && contains(lower(d.CurrentMode),'_ac')
                valid = isfield(d,'ACTimingMode') && ...
                    strcmp(d.ACTimingMode,'physical-time-fixed-rounds') && ...
                    isfield(d,'ACCommunicationInterval') && ...
                    abs(d.ACCommunicationInterval-0.01)<1e-12 && ...
                    isfield(d,'ACChecksPerPhysicalStep') && ...
                    d.ACChecksPerPhysicalStep==1;
            end
        case 'tp'
            valid = strcmp(d.TPQueryPolicy,'own-agent-query') && ...
                isequal(d.TPConsensusStateSize,[4,6]) && ...
                d.TPDACFixedRounds==10 && d.TPACFixedRounds==10 && ...
                abs(d.TPACConsensusStep-0.2)<1e-12 && ...
                isfield(d,'TPImplementationVersion') && ...
                d.TPImplementationVersion==3 && ...
                strcmp(d.TPConsensusTiming, ...
                'dac-continuous-ac-inner-fixed');
        case 'aggregation'
            valid = aggregation_parameters_are_current(d);
        case 'offline'
            valid = all(d.OfflineDataQuantity_set==cfg.offline_per_agent);
    end
    if valid && ismember(kind,{'ip','tp'})
        valid = aggregation_parameters_are_current(d);
    end
catch
    valid = false;
end

function valid = aggregation_parameters_are_current(d)
expected = control_aggregation_parameters();
if isfield(d,'AggregationParameters')
    actual = d.AggregationParameters;
elseif isfield(d,'aggregation_cfg')
    actual = d.aggregation_cfg;
else
    valid = false;
    return;
end
fields = {'posterior_var_floor','precision_floor','gpoe_beta_max', ...
    'rbcm_beta_max','bcm_prior_scale'};
valid = true;
for i = 1:numel(fields)
    name = fields{i};
    valid = valid && isfield(actual,name) && ...
        abs(actual.(name)-expected.(name)) < 1e-12;
end
end
end

function value = tracking_metric(result_file)
if ~isfile(result_file), error('Missing result: %s',result_file); end
d = load(result_file,'t_set','vartheta_all_set','TrackingError_vector');
if isfield(d,'vartheta_all_set') && ~isempty(d.vartheta_all_set)
    curve = sqrt(sum(d.vartheta_all_set.^2,1));
else
    curve = d.TrackingError_vector;
end
mask = d.t_set>=10;
value = max(curve(mask),[],'omitnan');
end

function value = communication_metric(result_file,kind)
if ~isfile(result_file), error('Missing result: %s',result_file); end
switch kind
    case 'ip_dac'
        d = load(result_file,'dac_broadcast_event_count_per_agent');
        value = mean(d.dac_broadcast_event_count_per_agent,'omitnan');
    case 'ip_ac'
        d = load(result_file,'ac_event_count_per_agent');
        value = mean(d.ac_event_count_per_agent,'omitnan');
    case 'tp'
        d = load(result_file,'tp_broadcast_count_per_agent');
        value = mean(d.tp_broadcast_count_per_agent,'omitnan');
    case 'nbr'
        d = load(result_file,'t_set');
        value = numel(d.t_set)-1;
    otherwise
        error('Unknown communication metric kind: %s',kind);
end
end

function write_complete_markdown(result_folder,tracking_mean,tracking_std, ...
    communication_mean,communication_std,methods,mode_names,seeds,cfg)
output_file = fullfile(result_folder,'complete_tracking_communication.md');
fid = fopen(output_file,'w');
if fid < 0, error('Cannot open Markdown output: %s',output_file); end
cleanup_file = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid,'# Complete control comparison\n\n');
fprintf(fid,['Seeds: %s. Formation: %d. T = %.0f s. Primary metric: ' ...
    '`max_{t >= 10 s} ||vartheta(t)||_2`. Values are Monte Carlo ' ...
    'mean +/- sample standard deviation.\n\n'], ...
    mat2str(seeds),cfg.use_formation,cfg.t_end);
fprintf(fid,'## Tracking error\n\n| Method |');
for mode_i = 1:numel(mode_names)
    fprintf(fid,' %s |',strrep(mode_names{mode_i},'_','-'));
end
fprintf(fid,'\n|---|');
for mode_i = 1:numel(mode_names), fprintf(fid,'---:|'); end
fprintf(fid,'\n');
for method_i = 1:numel(methods)
    fprintf(fid,'| %s |',upper(methods{method_i}));
    for mode_i = 1:numel(mode_names)
        fprintf(fid,' %s |',format_tracking( ...
            tracking_mean{method_i,mode_i},tracking_std{method_i,mode_i}));
    end
    fprintf(fid,'\n');
end

fprintf(fid,'\n## Communication statistics\n\n');
fprintf(fid,['Agent-level packaged broadcast instances per agent over ' ...
    'the complete trajectory.\n\n| Method |']);
for mode_i = 1:numel(mode_names)
    fprintf(fid,' %s |',strrep(mode_names{mode_i},'_','-'));
end
fprintf(fid,'\n|---|');
for mode_i = 1:numel(mode_names), fprintf(fid,'---:|'); end
fprintf(fid,'\n');
for method_i = 1:numel(methods)
    fprintf(fid,'| %s |',upper(methods{method_i}));
    for mode_i = 1:numel(mode_names)
        mean_value = communication_mean{method_i,mode_i};
        std_value = communication_std{method_i,mode_i};
        if isfinite(mean_value)
            fprintf(fid,' %.2f +/- %.2f |',mean_value,std_value);
        else
            fprintf(fid,' -- |');
        end
    end
    fprintf(fid,'\n');
end
fprintf(fid,['\nCEN, Offline and Exact do not have a distributed ' ...
    'agent-level broadcast count.\n']);
end

function text = format_tracking(mean_value,std_value)
if ~isfinite(mean_value)
    text = '--';
elseif std_value > 0 && std_value < 1e-5
    text = sprintf('%.5f +/- %.2e',mean_value,std_value);
else
    text = sprintf('%.5f +/- %.5f',mean_value,std_value);
end
end

function cfg = final_configuration(repo_root,use_formation)
cfg.t_end = 20;
cfg.unknown_scale = 0.4;
cfg.disturbance_scale = 0.1;
cfg.offline_per_agent = 300;
cfg.use_formation = logical(use_formation);
cfg.inducing_file = fullfile(repo_root,'result', ...
    'inducing_points_position100_velocity6_reference_seed42.mat');
if ~isfile(cfg.inducing_file)
    error('Missing structured inducing-point file: %s',cfg.inducing_file);
end
end

function configure_environment(cfg)
setenv('CONTROL_SIM_END_TIME',num2str(cfg.t_end));
setenv('CONTROL_ONLINE_AGENT_POLICY','all_agents');
setenv('CONTROL_IP_NUM_INDUCING_POINTS','600');
setenv('CONTROL_IP_INDUCING_POINT_FILE',cfg.inducing_file);
setenv('CONTROL_IP_POSITION_LENGTH_SCALE','0.5');
setenv('CONTROL_IP_VELOCITY_LENGTH_SCALE','1.5');
setenv('CONTROL_IP_RECONSTRUCTION_SIGMA_N','0.000001');
setenv('CONTROL_IP_PROJECTION_DIAGNOSTICS','0');
setenv('CONTROL_DAC_STEPS_PER_PHYSICAL_STEP','1');
setenv('CONTROL_DAC_TRIGGER_EPSILON','0.005');
setenv('CONTROL_AC_FIXED_ITERATIONS','10');
setenv('CONTROL_AC_CONSENSUS_STEP','0.2');
setenv('CONTROL_AC_KAPPA_P','20');
setenv('CONTROL_AC_PERIODIC_SIGMA','0.3');
setenv('CONTROL_AC_FORCE_FULL_BROADCAST','0');
setenv('CONTROL_AC_TRIGGER_DIAGNOSTICS','0');
setenv('CONTROL_AGGREGATION_VARIANCE_FLOOR','0.001');
setenv('CONTROL_AGGREGATION_PRECISION_FLOOR','0.0001');
setenv('CONTROL_GPOE_BETA_MAX','10');
setenv('CONTROL_RBCM_BETA_MAX','0.25');
setenv('CONTROL_BCM_PRIOR_SCALE','0.25');
setenv('CONTROL_TP_DAC_FIXED_ROUNDS','10');
setenv('CONTROL_TP_QUERY_UPDATE_INTERVAL','0.1');
setenv('CONTROL_TP_AC_FIXED_ROUNDS','10');
setenv('CONTROL_TP_AC_CONSENSUS_STEP','0.2');
end
