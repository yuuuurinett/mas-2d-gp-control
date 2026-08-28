%% RUN_CONTROL_MC10
% Monte Carlo runner and Word-friendly IP-only summary exporter.
%
% MC trials change only the agent initial-position seed.  This time-limited
% summary compares IP-DAC/IP-AC online and offline control modes.

clc;
script_folder = fileparts(mfilename('fullpath'));
if isempty(script_folder)
    script_folder = fileparts(which('run_control_mc10'));
end
project_root = fileparts(script_folder);
addpath(genpath(project_root));

%% Configuration
use_formation = true;
form_tag_set = {'noformation', 'formation'};
form_tag = form_tag_set{double(use_formation) + 1};

mc_count = 10;
initial_seed_set = 1001:(1000 + mc_count);
SimulationHorizonTag = 'T10_M800_th035';
tag_override = strtrim(getenv('CONTROL_MC10_RESULT_TAG'));
if ~isempty(tag_override)
    SimulationHorizonTag = tag_override;
end
StableStartTime = 3.0;

% Reproducible control-task configuration selected by the single-seed
% threshold sensitivity check.  Set these here so MC10 does not depend on
% stale environment variables left in the MATLAB session.
setenv('CONTROL_SIM_END_TIME', '');
setenv('CONTROL_IP_NUM_INDUCING_POINTS', '800');
setenv('CONTROL_ONLINE_TRIGGER_ERROR_TOL', '0.05');
setenv('CONTROL_CONSENSUS_INPUT_TRIGGER_TOL', '0.35');
setenv('CONTROL_DAC_CONSENSUS_INPUT_TRIGGER_TOL', '');
setenv('CONTROL_AC_CONSENSUS_INPUT_TRIGGER_TOL', '');
setenv('CONTROL_PLOT_INNER_CONSENSUS_TRIGGER', '');
if isempty(getenv('CONTROL_RBCM_BETA_MAX'))
    setenv('CONTROL_RBCM_BETA_MAX', '0.75');
end

UnknownScaleOverride = [];
DisturbanceScaleOverride = [];
OfflineDataQuantity = 350; % per agent, as in the shared manipulator code
run_missing_trials = strcmp(getenv('RUN_CONTROL_MC10_RUN_MISSING'), '1');

agg_methods = {'poe','gpoe','moe','bcm','rbcm'};
method_filter_raw = strtrim(getenv('CONTROL_MC10_METHODS'));
if ~isempty(method_filter_raw)
    requested_methods = regexp(lower(method_filter_raw), '[,\s;]+', 'split');
    requested_methods = requested_methods(~cellfun('isempty', requested_methods));
    invalid_methods = setdiff(requested_methods, agg_methods);
    if ~isempty(invalid_methods)
        error('Unknown CONTROL_MC10_METHODS entry: %s', strjoin(invalid_methods, ', '));
    end
    agg_methods = requested_methods;
end

run_plan = struct('key', {}, 'runner', {}, 'mode', {});
for method_idx = 1:numel(agg_methods)
    method = agg_methods{method_idx};
    ac_method = [method, '_ac'];

    run_plan(end+1) = make_plan(method,                  'ip',  method); %#ok<SAGROW>
    run_plan(end+1) = make_plan(ac_method,               'ip',  ac_method); %#ok<SAGROW>
    run_plan(end+1) = make_plan([method, '_offline'],    'ip',  [method, '_offline']); %#ok<SAGROW>
    run_plan(end+1) = make_plan([ac_method, '_offline'], 'ip',  [ac_method, '_offline']); %#ok<SAGROW>
end

result_folder_name = ['Control_MC10_', SimulationHorizonTag];
project_result_folder = fullfile(project_root, 'result', result_folder_name);
legacy_result_folder = fullfile(script_folder, 'result', result_folder_name);
if isfolder(legacy_result_folder) && ...
        ~isempty(dir(fullfile(legacy_result_folder, '*.mat'))) && ...
        isempty(dir(fullfile(project_result_folder, '*.mat')))
    result_folder = legacy_result_folder;
else
    result_folder = project_result_folder;
end
table_folder = fullfile(result_folder, 'Tables');
if ~isfolder(result_folder)
    mkdir(result_folder);
end
if ~isfolder(table_folder)
    mkdir(table_folder);
end

%% Run MC trials
mode_count = numel(run_plan);
final_error_matrix = nan(mode_count, mc_count);
stable_worst_error_matrix = nan(mode_count, mc_count);
stable_agent_avg_error_matrix = nan(mode_count, mc_count);
per_agent_stable_max_matrix = nan(mode_count, mc_count, 6);
trigger_matrix = nan(mode_count, mc_count);
online_update_matrix = nan(mode_count, mc_count);
elapsed_matrix = nan(mode_count, mc_count);

for mode_idx = 1:mode_count
    plan_item = run_plan(mode_idx);
    for mc_idx = 1:mc_count
        initial_seed = initial_seed_set(mc_idx);
        save_name = sprintf('%s_%s_mc%02d_seed%d', ...
            plan_item.key, form_tag, mc_idx, initial_seed);
        save_file = fullfile(result_folder, [save_name, '.mat']);

        if isfile(save_file)
            fprintf('[SKIP] %s\n', save_name);
        elseif ~run_missing_trials
            fprintf('[MISSING] %s\n', save_name);
        else
            fprintf('[RUN] key=%s, runner=%s, mode=%s, MC=%d/%d, seed=%d\n', ...
                plan_item.key, plan_item.runner, plan_item.mode, ...
                mc_idx, mc_count, initial_seed);
            try
                run_one_control_trial(plan_item, result_folder, save_name, ...
                    use_formation, UnknownScaleOverride, DisturbanceScaleOverride, ...
                    OfflineDataQuantity, initial_seed);
            catch ME
                error_file = fullfile(result_folder, [save_name, '_ERROR.mat']);
                error_message = getReport(ME, 'extended', 'hyperlinks', 'off');
                save(error_file, 'plan_item', 'mc_idx', 'initial_seed', ...
                    'error_message', 'ME');
                warning('MC trial failed and was recorded: %s\n%s', ...
                    error_file, ME.message);
            end
        end

        error_file = fullfile(result_folder, [save_name, '_ERROR.mat']);
        stats = load_control_mc_stats(save_file, error_file, StableStartTime);
        final_error_matrix(mode_idx, mc_idx) = stats.final_error;
        stable_worst_error_matrix(mode_idx, mc_idx) = stats.stable_worst_error;
        stable_agent_avg_error_matrix(mode_idx, mc_idx) = stats.stable_agent_avg_error;
        agent_count_stats = min(numel(stats.per_agent_stable_max), ...
            size(per_agent_stable_max_matrix, 3));
        if agent_count_stats > 0
            per_agent_stable_max_matrix(mode_idx, mc_idx, 1:agent_count_stats) = ...
                stats.per_agent_stable_max(1:agent_count_stats);
        end
        trigger_matrix(mode_idx, mc_idx) = stats.trigger_avg;
        online_update_matrix(mode_idx, mc_idx) = stats.online_update_avg;
        elapsed_matrix(mode_idx, mc_idx) = stats.elapsed_time;
    end
end

%% Export summaries
mode_key_list = string({run_plan.key}');
summary = table(mode_key_list, ...
    mean(stable_worst_error_matrix, 2, 'omitnan'), ...
    std(stable_worst_error_matrix, 0, 2, 'omitnan'), ...
    mean(stable_agent_avg_error_matrix, 2, 'omitnan'), ...
    std(stable_agent_avg_error_matrix, 0, 2, 'omitnan'), ...
    mean(final_error_matrix, 2, 'omitnan'), ...
    std(final_error_matrix, 0, 2, 'omitnan'), ...
    sum(isnan(stable_worst_error_matrix), 2), ...
    mean(trigger_matrix, 2, 'omitnan'), ...
    std(trigger_matrix, 0, 2, 'omitnan'), ...
    mean(online_update_matrix, 2, 'omitnan'), ...
    mean(elapsed_matrix, 2, 'omitnan'), ...
    'VariableNames', {'Mode', ...
    'StableWorstErrorMean', 'StableWorstErrorStd', ...
    'StableAgentAvgErrorMean', 'StableAgentAvgErrorStd', ...
    'FinalErrorMean', 'FinalErrorStd', ...
    'FailedTrials', 'TriggerMean', 'TriggerStd', ...
    'OnlineUpdatesMean', 'ElapsedTimeMean'});

mat_file = fullfile(result_folder, sprintf('control_mc10_summary_%s.mat', form_tag));
csv_file = fullfile(table_folder, sprintf('control_mc10_summary_%s.csv', form_tag));
html_file = fullfile(table_folder, sprintf('control_mc10_summary_%s.html', form_tag));

save(mat_file, 'summary', 'run_plan', 'initial_seed_set', ...
    'SimulationHorizonTag', 'StableStartTime', 'final_error_matrix', ...
    'stable_worst_error_matrix', 'stable_agent_avg_error_matrix', ...
    'per_agent_stable_max_matrix', 'trigger_matrix', ...
    'online_update_matrix', 'elapsed_matrix');
writetable(summary, csv_file);
write_summary_html(html_file, summary, form_tag, mc_count, StableStartTime);

fprintf('\nMC10 summary written to:\n%s\n%s\n%s\n', mat_file, csv_file, html_file);

%% Local helper functions
function item = make_plan(key, runner, mode)
item = struct('key', key, 'runner', runner, 'mode', mode);
end

function run_one_control_trial(plan_item, result_folder, save_name, ...
    use_formation, unknown_scale, disturbance_scale, offline_data_quantity, initial_seed)
is_offline = contains(plan_item.mode, '_offline');
switch plan_item.runner
    case 'ip'
        if is_offline
            run_simulation_inducing_point(plan_item.mode, result_folder, save_name, ...
                use_formation, unknown_scale, disturbance_scale, ...
                offline_data_quantity, initial_seed);
        else
            run_simulation_inducing_point(plan_item.mode, result_folder, save_name, ...
                use_formation, unknown_scale, disturbance_scale, [], initial_seed);
        end
    case 'tp'
        if is_offline
            run_simulation_test_point(plan_item.mode, result_folder, save_name, ...
                use_formation, initial_seed, offline_data_quantity);
        else
            run_simulation_test_point(plan_item.mode, result_folder, save_name, ...
                use_formation, initial_seed);
        end
    case 'cen'
        if is_offline
            run_simulation_cen(plan_item.mode, result_folder, save_name, ...
                use_formation, initial_seed, offline_data_quantity);
        else
            run_simulation_cen(plan_item.mode, result_folder, save_name, ...
                use_formation, initial_seed);
        end
    case 'nbr'
        if is_offline
            run_simulation_nbr(plan_item.mode, result_folder, save_name, ...
                use_formation, initial_seed, offline_data_quantity);
        else
            run_simulation_nbr(plan_item.mode, result_folder, save_name, ...
                use_formation, initial_seed);
        end
    otherwise
        error('Unknown runner: %s', plan_item.runner);
end
end

function stats = load_control_mc_stats(save_file, error_file, stable_start_time)
stats = struct('final_error', NaN, 'stable_worst_error', NaN, ...
    'stable_agent_avg_error', NaN, 'per_agent_stable_max', NaN(1, 0), ...
    'trigger_avg', NaN, ...
    'online_update_avg', NaN, 'elapsed_time', NaN);
if nargin >= 2 && isfile(error_file) && ~isfile(save_file)
    return;
end
if ~isfile(save_file)
    return;
end

d = load(save_file);
if isfield(d, 'vartheta_all_set')
    stats.final_error = final_agent_average_tracking_error(d);
    [stats.stable_worst_error, stats.stable_agent_avg_error, ...
        stats.per_agent_stable_max] = stable_tracking_error_metrics( ...
        d, stable_start_time);
elseif isfield(d, 'TrackingError_vector')
    stats.final_error = d.TrackingError_vector(end);
    if isfield(d, 't_set')
        stable_mask = d.t_set >= stable_start_time;
        stable_error = d.TrackingError_vector(stable_mask);
        if any(isfinite(stable_error))
            stats.stable_worst_error = max(stable_error, [], 'omitnan');
            stats.stable_agent_avg_error = stats.stable_worst_error;
        end
    end
end
if isfield(d, 'online_update_count')
    stats.online_update_avg = mean(d.online_update_count(:), 'omitnan');
elseif isfield(d, 'online_trigger_count')
    stats.online_update_avg = mean(d.online_trigger_count(:), 'omitnan');
end
if isfield(d, 'elapsed_time')
    stats.elapsed_time = d.elapsed_time;
end

mode_name = "";
if isfield(d, 'CurrentMode')
    mode_name = lower(string(d.CurrentMode));
end

if contains(mode_name, "_ac") && isfield(d, 'ac_event_count_per_agent')
    stats.trigger_avg = mean(d.ac_event_count_per_agent(:), 'omitnan');
elseif isfield(d, 'dac_event_count_per_agent')
    stats.trigger_avg = mean(d.dac_event_count_per_agent(:), 'omitnan');
elseif isfield(d, 'dac_broadcasts_per_agent')
    stats.trigger_avg = d.dac_broadcasts_per_agent;
elseif isfield(d, 'ac_broadcasts_per_agent')
    stats.trigger_avg = d.ac_broadcasts_per_agent;
elseif isfield(d, 'online_trigger_count')
    stats.trigger_avg = mean(d.online_trigger_count(:), 'omitnan');
end
end

function [system_worst_error, agent_avg_error, per_agent_stable_max] = ...
    stable_tracking_error_metrics(d, stable_start_time)
system_worst_error = NaN;
agent_avg_error = NaN;
per_agent_stable_max = NaN(1, 0);

agent_error = compute_agent_tracking_error(d);
if isempty(agent_error) || ~isfield(d, 't_set')
    return;
end

t = d.t_set(:).';
sample_count = min(numel(t), size(agent_error, 2));
t = t(1:sample_count);
agent_error = agent_error(:,1:sample_count);
stable_mask = t >= stable_start_time;
if ~any(stable_mask)
    return;
end

stable_error = agent_error(:,stable_mask);
per_agent_stable_max = max(stable_error, [], 2, 'omitnan').';
finite_agent_mask = isfinite(per_agent_stable_max);
if ~any(finite_agent_mask)
    return;
end

system_worst_error = max(per_agent_stable_max(finite_agent_mask));
agent_avg_error = mean(per_agent_stable_max(finite_agent_mask), 'omitnan');
end

function agent_error = compute_agent_tracking_error(d)
agent_error = [];
if ~isfield(d, 'vartheta_all_set')
    return;
end

agent_quantity = infer_agent_quantity(d);
state_count = size(d.vartheta_all_set, 1);
if agent_quantity <= 0 || mod(state_count, agent_quantity) ~= 0
    return;
end

state_dim = state_count / agent_quantity;
sample_count = size(d.vartheta_all_set, 2);
agent_error = nan(agent_quantity, sample_count);
for sample_idx = 1:sample_count
    vartheta_sample = reshape(d.vartheta_all_set(:,sample_idx), ...
        state_dim, agent_quantity);
    agent_error(:,sample_idx) = vecnorm(vartheta_sample, 2, 1).';
end
end

function final_error = final_agent_average_tracking_error(d)
final_error = NaN;
vartheta_final = d.vartheta_all_set(:,end);
if any(~isfinite(vartheta_final))
    return;
end

agent_quantity = infer_agent_quantity(d);
if agent_quantity <= 0 || mod(numel(vartheta_final), agent_quantity) ~= 0
    return;
end

state_dim = numel(vartheta_final) / agent_quantity;
vartheta_by_agent = reshape(vartheta_final, state_dim, agent_quantity);
agent_error = vecnorm(vartheta_by_agent, 2, 1);
final_error = mean(agent_error, 'omitnan');
end

function agent_quantity = infer_agent_quantity(d)
if isfield(d, 'online_update_count') && ~isempty(d.online_update_count)
    agent_quantity = numel(d.online_update_count);
elseif isfield(d, 'dac_event_count_per_agent') && ~isempty(d.dac_event_count_per_agent)
    agent_quantity = numel(d.dac_event_count_per_agent);
elseif isfield(d, 'ac_event_count_per_agent') && ~isempty(d.ac_event_count_per_agent)
    agent_quantity = numel(d.ac_event_count_per_agent);
else
    agent_quantity = 6;
end
end

function write_summary_html(html_file, summary, form_tag, mc_count, stable_start_time)
fid = fopen(html_file, 'w');
if fid < 0
    error('Cannot write HTML summary: %s', html_file);
end
cleanup_obj = onCleanup(@() fclose(fid));

fprintf(fid, '<!DOCTYPE html>\n<html><head><meta charset="UTF-8">\n');
fprintf(fid, '<style>\n');
fprintf(fid, 'body{font-family:"Times New Roman",serif;margin:24px;color:#111;}\n');
fprintf(fid, 'h3{color:#c00000;margin:20px 0 6px 0;font-size:12pt;}\n');
fprintf(fid, 'table{border-collapse:collapse;font-size:11pt;margin-bottom:22px;}\n');
fprintf(fid, 'th,td{border:1px solid #555;padding:6px 9px;text-align:center;}\n');
fprintf(fid, 'th{background:#bfe3f1;font-weight:bold;}\n');
fprintf(fid, 'td.mode{font-weight:bold;background:#f7f7f7;}\n');
fprintf(fid, '.small{font-size:10pt;margin-top:-12px;}\n');
fprintf(fid, '</style></head><body>\n');

methods = {'poe','gpoe','moe','bcm','rbcm'};
method_labels = {'POE','GPOE','MOE','BCM','RBCM'};

fprintf(fid, '<h3>Stable-window worst tracking error - MC%d mean [%s]</h3>\n', mc_count, form_tag);
fprintf(fid, '<table>\n');
fprintf(fid, '<tr><th>Method</th><th>IP-DAC online</th><th>IP-AC online</th><th>IP-DAC offline</th><th>IP-AC offline</th></tr>\n');
for method_idx = 1:numel(methods)
    method = methods{method_idx};
    fprintf(fid, '<tr>');
    fprintf(fid, '<td class="mode">%s</td>', method_labels{method_idx});
    fprintf(fid, '<td>%s</td>', format_error_cell(summary, method));
    fprintf(fid, '<td>%s</td>', format_error_cell(summary, [method, '_ac']));
    fprintf(fid, '<td>%s</td>', format_error_cell(summary, [method, '_offline']));
    fprintf(fid, '<td>%s</td>', format_error_cell(summary, [method, '_ac_offline']));
    fprintf(fid, '</tr>\n');
end
fprintf(fid, '</table>\n');

fprintf(fid, '<h3>Communication statistics [%s]</h3>\n', form_tag);
fprintf(fid, '<table>\n');
fprintf(fid, '<tr><th>Method</th><th>IP-DAC online<br/>comm. / agent</th><th>IP-AC online<br/>comm. / agent</th><th>IP-DAC offline<br/>comm. / agent</th><th>IP-AC offline<br/>comm. / agent</th></tr>\n');
for method_idx = 1:numel(methods)
    method = methods{method_idx};
    fprintf(fid, '<tr>');
    fprintf(fid, '<td class="mode">%s</td>', method_labels{method_idx});
    fprintf(fid, '<td>%s</td>', format_trigger_cell(summary, method));
    fprintf(fid, '<td>%s</td>', format_trigger_cell(summary, [method, '_ac']));
    fprintf(fid, '<td>%s</td>', format_trigger_cell(summary, [method, '_offline']));
    fprintf(fid, '<td>%s</td>', format_trigger_cell(summary, [method, '_ac_offline']));
    fprintf(fid, '</tr>\n');
end
fprintf(fid, '</table>\n');

fprintf(fid, '<p class="small">Each tracking-error entry is the MC%d mean of the MAS worst-case stable-window error: first compute each agent''s e_i(t)=||vartheta_i(t)||, take max over t &ge; %.1f s for each agent, then take max over agents. ', mc_count, stable_start_time);
fprintf(fid, 'The saved MAT file also includes the per-agent stable maxima and the agent-averaged stable maxima. ');
fprintf(fid, 'Only IP-DAC/IP-AC results are included in this time-limited control summary. ');
fprintf(fid, 'The without-GP no-compensation baseline is excluded from MC averaging and should be shown separately as a divergence case if needed. ');
fprintf(fid, 'Cells marked with "-" were not run or are unavailable; cells marked with "fail=n" contain non-finite or missing trials.</p>\n');
fprintf(fid, '</body></html>\n');
end

function txt = format_error_cell(summary, mode_name)
row_idx = find(summary.Mode == string(mode_name), 1);
if isempty(row_idx) || ~isfinite(summary.StableWorstErrorMean(row_idx))
    txt = '-';
    return;
end
txt = sprintf('%.4f', summary.StableWorstErrorMean(row_idx));
if summary.FailedTrials(row_idx) > 0
    txt = sprintf('%s (fail=%d)', txt, summary.FailedTrials(row_idx));
end
end

function txt = format_trigger_cell(summary, mode_name)
row_idx = find(summary.Mode == string(mode_name), 1);
if isempty(row_idx) || ~isfinite(summary.TriggerMean(row_idx))
    txt = '-';
    return;
end
txt = sprintf('%.2f', summary.TriggerMean(row_idx));
if summary.FailedTrials(row_idx) > 0
    txt = sprintf('%s (fail=%d)', txt, summary.FailedTrials(row_idx));
end
end
