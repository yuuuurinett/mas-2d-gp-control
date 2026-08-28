function [MeanTable,StdTable,LongTable] = ...
    run_mc10_all_aggregation_methods( ...
    do_run,seed_subset,do_summarize,num_inducing_points,ip_only,mc_count)
%RUN_MC10_ALL_AGGREGATION_METHODS Configurable-seed control comparison.
% The historical function name is retained for compatibility.  Use the
% sixth input to request a different Monte Carlo count, or call the MC50
% wrapper run_mc50_all_aggregation_methods.
% Metric: max_{t>=10 s} ||vartheta_all(t)||_2, matching the advisor's
% manipulator-code convention. Online modes collect one local sample per
% agent every 0.1 s; Local is offline with 350 samples per agent, matching
% the advisor's manipulator shared code.

if nargin < 1 || isempty(do_run), do_run = true; end
if nargin < 6 || isempty(mc_count), mc_count = 10; end
mc_count = max(1,round(mc_count));
if mc_count == 10
    all_seed_values = 42:51;
else
    % New larger batches use a self-contained 1:N seed convention.
    all_seed_values = 1:mc_count;
end
if nargin < 2 || isempty(seed_subset)
    seed_values_to_run = all_seed_values;
else
    seed_values_to_run = seed_subset(:).';
end
if nargin < 3 || isempty(do_summarize), do_summarize = true; end
if nargin < 4 || isempty(num_inducing_points), num_inducing_points = 600; end
num_inducing_points = max(1,round(num_inducing_points));
if nargin < 5 || isempty(ip_only)
    ip_only = num_inducing_points ~= 400;
end
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result', ...
    sprintf('mc%d_all_aggregation_m%d_T30_per_agent', ...
    mc_count,num_inducing_points));
if ~isfolder(result_folder), mkdir(result_folder); end
shared_result_folder = fullfile(repo_root,'result', ...
    sprintf('mc%d_all_aggregation_m400_T30_per_agent',mc_count));
if ip_only && ~isfolder(shared_result_folder)
    error('Shared M400 baseline folder is missing: %s',shared_result_folder);
end
legacy_result_folder = '';
% Do not copy MC10 data into larger batches.  The IP sampling, Offline
% allocation and dynamics settings evolved after those files were made;
% mixing them would invalidate the Monte Carlo statistics.

aggregation_methods = {'poe','gpoe','moe','bcm','rbcm'};
mode_names = {'IP_DAC','IP_AC','TP_DAC','TP_AC','CEN','NBR', ...
    'Offline','Exact'};
seed_values = all_seed_values;
offline_allocation = 350;

if do_run
    for seed_idx = 1:numel(seed_values_to_run)
        run_one_seed(seed_values_to_run(seed_idx),aggregation_methods, ...
            offline_allocation,result_folder,num_inducing_points,ip_only, ...
            legacy_result_folder);
    end
end

if ~do_summarize
    MeanTable = table(); StdTable = table(); LongTable = table();
    return;
end

method_count = numel(aggregation_methods);
mode_count = numel(mode_names);
seed_count = numel(seed_values);
metric = nan(method_count,mode_count,seed_count);

for seed_idx = 1:seed_count
    seed = seed_values(seed_idx);
    seed_folder = fullfile(result_folder,sprintf('seed_%03d',seed));
    common_seed_folder = seed_folder;
    if ip_only
        common_seed_folder = fullfile(shared_result_folder, ...
            sprintf('seed_%03d',seed));
    end
    local_value = tracking_metric(fullfile( ...
        common_seed_folder,'local_offline_peragent350_seeded.mat'));
    exact_value = tracking_metric(fullfile(common_seed_folder,'exact.mat'));
    for method_idx = 1:method_count
        method = aggregation_methods{method_idx};
        metric(method_idx,1,seed_idx) = tracking_metric( ...
            fullfile(seed_folder,[method '_ip_dac.mat']));
        metric(method_idx,2,seed_idx) = tracking_metric( ...
            fullfile(seed_folder,[method '_ip_ac.mat']));
        metric(method_idx,3,seed_idx) = tracking_metric( ...
            fullfile(common_seed_folder,[method '_tp_dac_fixedR10.mat']));
        metric(method_idx,4,seed_idx) = tracking_metric( ...
            fullfile(common_seed_folder,[method '_tp_ac_fixedR10.mat']));
        metric(method_idx,5,seed_idx) = tracking_metric( ...
            fullfile(common_seed_folder,[method '_cen.mat']));
        metric(method_idx,6,seed_idx) = tracking_metric( ...
            fullfile(common_seed_folder,[method '_nbr.mat']));
        metric(method_idx,7,seed_idx) = local_value;
        metric(method_idx,8,seed_idx) = exact_value;
    end
end

mean_metric = mean(metric,3,'omitnan');
std_metric = std(metric,0,3,'omitnan');
aggregation_row_names = cellstr(upper(string(aggregation_methods)));
MeanTable = array2table(mean_metric,'VariableNames',mode_names, ...
    'RowNames',aggregation_row_names);
StdTable = array2table(std_metric,'VariableNames',mode_names, ...
    'RowNames',aggregation_row_names);
formatted_metric = strings(method_count,mode_count);
for method_idx = 1:method_count
    for mode_idx = 1:mode_count
        if std_metric(method_idx,mode_idx) > 0 && ...
                std_metric(method_idx,mode_idx) < 1e-5
            formatted_metric(method_idx,mode_idx) = sprintf( ...
                '%.5f +/- %.2e',mean_metric(method_idx,mode_idx), ...
                std_metric(method_idx,mode_idx));
        else
            formatted_metric(method_idx,mode_idx) = sprintf( ...
                '%.5f +/- %.5f',mean_metric(method_idx,mode_idx), ...
                std_metric(method_idx,mode_idx));
        end
    end
end
MeanStdTable = array2table(formatted_metric,'VariableNames',mode_names, ...
    'RowNames',aggregation_row_names);
writetable(MeanTable,fullfile(result_folder, ...
    'tracking_max_after10s_mean.csv'),'WriteRowNames',true);
writetable(StdTable,fullfile(result_folder, ...
    'tracking_max_after10s_std.csv'),'WriteRowNames',true);
writetable(MeanStdTable,fullfile(result_folder, ...
    'tracking_max_after10s_mean_std.csv'),'WriteRowNames',true);

row_count = method_count*mode_count*seed_count;
Aggregation = strings(row_count,1);
Mode = strings(row_count,1);
Seed = zeros(row_count,1);
MaxTrackingErrorAfter10s = nan(row_count,1);
row = 0;
for method_idx = 1:method_count
    for mode_idx = 1:mode_count
        for seed_idx = 1:seed_count
            row = row+1;
            Aggregation(row) = upper(aggregation_methods{method_idx});
            Mode(row) = mode_names{mode_idx};
            Seed(row) = seed_values(seed_idx);
            MaxTrackingErrorAfter10s(row) = ...
                metric(method_idx,mode_idx,seed_idx);
        end
    end
end
LongTable = table(Aggregation,Mode,Seed,MaxTrackingErrorAfter10s);
writetable(LongTable,fullfile(result_folder, ...
    'tracking_max_after10s_all_seeds.csv'));
save(fullfile(result_folder,sprintf('mc%d_summary.mat',mc_count)), ...
    'MeanTable','StdTable', ...
    'MeanStdTable', ...
    'LongTable','metric','aggregation_methods','mode_names','seed_values', ...
    'offline_allocation','num_inducing_points','ip_only', ...
    'shared_result_folder');
disp('Mean max tracking error for t >= 10 s:');
disp(MeanTable);
fprintf('Standard deviation across %d seeds:\n',mc_count);
disp(StdTable);
disp('Mean +/- standard deviation (5 decimals; scientific notation when needed):');
disp(MeanStdTable);
end

function run_one_seed(seed,aggregation_methods,offline_allocation, ...
    result_folder,num_inducing_points,ip_only,legacy_result_folder)
configure_worker_environment(num_inducing_points);
seed_folder = fullfile(result_folder,sprintf('seed_%03d',seed));
if ~isfolder(seed_folder), mkdir(seed_folder); end
copy_legacy_seed_results(legacy_result_folder,seed_folder,seed);

for method_idx = 1:numel(aggregation_methods)
    method = aggregation_methods{method_idx};
    run_quiet_if_missing(fullfile(seed_folder,[method '_ip_dac.mat']), ...
        @() run_simulation_inducing_point(method,seed_folder, ...
        [method '_ip_dac'],false,[],[],[],seed));
    run_quiet_if_missing(fullfile(seed_folder,[method '_ip_ac.mat']), ...
        @() run_simulation_inducing_point([method '_ac'],seed_folder, ...
        [method '_ip_ac'],false,[],[],[],seed));
    if ip_only
        continue;
    end
    run_quiet_if_missing(fullfile(seed_folder, ...
        [method '_tp_dac_fixedR10.mat']), ...
        @() run_simulation_test_point(method,seed_folder, ...
        [method '_tp_dac_fixedR10'],false,seed,0));
    run_quiet_if_missing(fullfile(seed_folder, ...
        [method '_tp_ac_fixedR10.mat']), ...
        @() run_simulation_test_point([method '_ac'],seed_folder, ...
        [method '_tp_ac_fixedR10'],false,seed,0));
    run_quiet_if_missing(fullfile(seed_folder,[method '_cen.mat']), ...
        @() run_simulation_cen(method,seed_folder,[method '_cen'], ...
        false,seed,0));
    run_quiet_if_missing(fullfile(seed_folder,[method '_nbr.mat']), ...
        @() run_simulation_nbr(method,seed_folder,[method '_nbr'], ...
        false,seed,0));
end

if ~ip_only
    run_quiet_if_missing(fullfile(seed_folder,'local_offline_peragent350_seeded.mat'), ...
        @() run_simulation_test_point('local_offline',seed_folder, ...
        'local_offline_peragent350_seeded',false,seed,offline_allocation));
    run_quiet_if_missing(fullfile(seed_folder,'exact.mat'), ...
        @() run_simulation_test_point('exact',seed_folder,'exact', ...
        false,seed,0));
end
if ip_only
    fprintf('Monte Carlo seed %d IP-DAC/IP-AC complete.\n',seed);
else
    fprintf('Monte Carlo seed %d complete.\n',seed);
end
end

function copy_legacy_seed_results(legacy_result_folder,seed_folder,seed)
% Reuse the completed first ten trials when extending MC10 to MC50.
if isempty(legacy_result_folder), return; end
legacy_seed_folder = fullfile(legacy_result_folder, ...
    sprintf('seed_%03d',seed));
if ~isfolder(legacy_seed_folder), return; end
legacy_files = dir(fullfile(legacy_seed_folder,'*.mat'));
for file_idx = 1:numel(legacy_files)
    source_file = fullfile(legacy_files(file_idx).folder, ...
        legacy_files(file_idx).name);
    target_file = fullfile(seed_folder,legacy_files(file_idx).name);
    if ~isfile(target_file), copyfile(source_file,target_file); end
end
end

function configure_worker_environment(num_inducing_points)
setenv('CONTROL_SIM_END_TIME','30');
setenv('CONTROL_ONLINE_AGENT_POLICY','all_agents');
setenv('CONTROL_IP_NUM_INDUCING_POINTS',num2str(num_inducing_points));
setenv('CONTROL_IP_INDUCING_POINT_FILE','');
setenv('CONTROL_IP_PROJECTION_DIAGNOSTICS','0');
setenv('CONTROL_DAC_STEPS_PER_PHYSICAL_STEP','1');
setenv('CONTROL_DAC_TRIGGER_EPSILON','0.005');
setenv('CONTROL_AC_FIXED_ITERATIONS','10');
setenv('CONTROL_AC_PERIODIC_SIGMA','0.3');
setenv('CONTROL_AC_TRIGGER_DIAGNOSTICS','0');
setenv('CONTROL_TP_DAC_STEPS_PER_PHYSICAL_STEP','1');
setenv('CONTROL_TP_DAC_FIXED_ROUNDS','10');
setenv('CONTROL_TP_DAC_INITIAL_ROUNDS','1');
setenv('CONTROL_TP_QUERY_UPDATE_INTERVAL','0.1');
setenv('CONTROL_TP_AC_ITERATION_POLICY','fixed');
setenv('CONTROL_TP_AC_FIXED_ROUNDS','10');
end

function run_quiet_if_missing(target,runner)
if isfile(target), return; end
evalc('runner();');
end

function value = tracking_metric(result_file)
if ~isfile(result_file), error('Missing result: %s',result_file); end
d = load(result_file,'t_set','vartheta_all_set','TrackingError_vector');
if isfield(d,'vartheta_all_set') && ~isempty(d.vartheta_all_set)
    tracking = sqrt(sum(d.vartheta_all_set.^2,1));
else
    tracking = d.TrackingError_vector(:).';
end
mask = d.t_set >= 10;
value = max(tracking(mask),[],'omitnan');
end
