function plot_time_vs_trigger_ip(method_list, form_tag, result_folder, output_folder)
%PLOT_TIME_VS_TRIGGER_IP Plot IP-DAC/IP-AC per-agent broadcast triggers.
%
% This figure is intended for advisor-facing control-task results.  It
% visualizes IP-DAC per-agent communication events on the physical simulation
% time axis, similar to the trigger-instance plot in the cooperative online
% learning paper and the shared manipulator examples.
%
% Usage:
%   plot_time_vs_trigger_ip
%   plot_time_vs_trigger_ip('poe')
%   plot_time_vs_trigger_ip({'poe','bcm','rbcm'}, 'formation')

if nargin < 1 || isempty(method_list)
    method_list = {'poe','gpoe','moe','bcm','rbcm'};
end
if ischar(method_list) || isstring(method_list)
    method_list = cellstr(method_list);
end
if nargin < 2 || isempty(form_tag)
    form_tag = 'formation';
end
if nargin < 3 || isempty(result_folder)
    script_folder = fileparts(mfilename('fullpath'));
    result_folder = fullfile(script_folder, 'result', 'Diagnostics');
end
if nargin < 4 || isempty(output_folder)
    script_folder = fileparts(mfilename('fullpath'));
    output_folder = fullfile(script_folder, 'result', 'Figures');
end
if ~isfolder(output_folder)
    mkdir(output_folder);
end
use_inner_consensus_trigger = strcmp(getenv('CONTROL_PLOT_INNER_CONSENSUS_TRIGGER'), '1');

for method_idx = 1:numel(method_list)
    method = lower(method_list{method_idx});
    dac_file = find_result_file(result_folder, method, form_tag);
    ac_method = [method, '_ac'];
    ac_file = find_result_file(result_folder, ac_method, form_tag);
    diagnostics_folder = fullfile(fileparts(mfilename('fullpath')), 'result', 'Diagnostics');

    if use_inner_consensus_trigger
        dac = load_inner_trigger_history(dac_file, 'dac_inner_trigger_instance_set');
    else
        dac = empty_trigger_history();
    end
    if isempty(dac.event_times)
        dac = load_trigger_history(dac_file, 'dac_broadcast_count_set');
    end
    if use_inner_consensus_trigger && isempty(dac.counts) && isempty(dac.event_times)
        dac = load_trigger_history_from_diagnostics(diagnostics_folder, method, form_tag, ...
            'dac_inner_trigger_instance_set', false, true);
    end
    if isempty(dac.counts) && isempty(dac.event_times)
        dac = load_trigger_history_from_diagnostics(diagnostics_folder, method, form_tag, ...
            'dac_broadcast_count_set', false, false);
    end

    if use_inner_consensus_trigger
        ac = load_inner_trigger_history(ac_file, 'ac_inner_trigger_instance_set');
    else
        ac = empty_trigger_history();
    end
    if isempty(ac.event_times)
        ac = load_trigger_history(ac_file, 'ac_physical_trigger_count_set');
    end
    if use_inner_consensus_trigger && isempty(ac.counts) && isempty(ac.event_times)
        ac = load_trigger_history_from_diagnostics(diagnostics_folder, ac_method, form_tag, ...
            'ac_inner_trigger_instance_set', true, true);
    end
    if isempty(ac.counts) && isempty(ac.event_times)
        ac = load_trigger_history_from_diagnostics(diagnostics_folder, ac_method, form_tag, ...
            'ac_physical_trigger_count_set', true, false);
    end

    if isempty(dac.counts) && isempty(dac.event_times) && ...
            isempty(ac.counts) && isempty(ac.event_times)
        warning('No IP-DAC/IP-AC physical-time trigger history found for method "%s".', method);
        continue;
    end

    fig = figure('Name', sprintf('Time-vs-trigger %s %s', upper(method), form_tag), ...
        'Color', 'w', 'Position', [100 100 980 620]);
    layout = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax_dac = nexttile(layout);
    plot_trigger_instances(ax_dac, dac, 'IP-DAC');
    title(ax_dac, 'IP-DAC', 'Interpreter', 'none');

    ax_ac = nexttile(layout);
    plot_trigger_instances(ax_ac, ac, 'IP-AC');
    title(ax_ac, 'IP-AC', 'Interpreter', 'none');

    if use_inner_consensus_trigger
        figure_title = sprintf('%s internal consensus broadcast trigger instances [%s]', ...
            upper(method), form_tag);
    else
        figure_title = sprintf('%s physical-time broadcast trigger instances [%s]', ...
            upper(method), form_tag);
    end
    title(layout, figure_title, ...
        'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold', ...
        'Interpreter', 'none');

    base_name = sprintf('TimeVsTrigger_IP_DAC_AC_%s_%s', upper(method), form_tag);
    png_file = fullfile(output_folder, [base_name, '.png']);
    if isfile(png_file)
        delete(png_file);
    end
    exportgraphics(fig, png_file, 'Resolution', 200);
    savefig(fig, fullfile(output_folder, [base_name, '.fig']));
end
end

function file_path = find_result_file(result_folder, method, form_tag)
file_path = fullfile(result_folder, sprintf('%s_%s.mat', method, form_tag));
if isfile(file_path)
    return;
end

pattern = sprintf('%s_%s_mc*_seed*.mat', method, form_tag);
candidates = dir(fullfile(result_folder, pattern));
diagnostic_pattern = sprintf('diag_%s*%s*.mat', method, form_tag);
diagnostic_candidates = [
    dir(fullfile(result_folder, diagnostic_pattern));
    dir(fullfile(result_folder, sprintf('diag_%s*.mat', method)))];
if ~contains(method, '_ac')
    diagnostic_candidates = diagnostic_candidates(~contains({diagnostic_candidates.name}, '_ac_'));
end
candidates = [candidates; diagnostic_candidates(:)];
if isempty(candidates)
    file_path = '';
    return;
end
[~, order] = sort([candidates.datenum], 'descend');
candidates = candidates(order);
file_path = fullfile(candidates(1).folder, candidates(1).name);
end

function history = load_trigger_history(file_path, field_name)
history = empty_trigger_history();
if ~isfile(file_path)
    return;
end

names = string({whos('-file', file_path).name});
if ~any(names == field_name) || ~any(names == "t_set")
    return;
end

data = load(file_path, field_name, 't_set');
counts = data.(field_name);
if isempty(counts)
    return;
end

sample_count = size(counts, 2);
time = data.t_set(:).';
history.counts = counts;
history.time = time(1:sample_count);
history.source = file_path;
end

function history = load_inner_trigger_history(file_path, field_name)
history = empty_trigger_history();
if ~isfile(file_path)
    return;
end

names = string({whos('-file', file_path).name});
if ~any(names == field_name) || ~any(names == "t_set")
    return;
end

load_names = {field_name, 't_set'};
if any(names == "IPOnlineUpdateInterval")
    load_names{end+1} = 'IPOnlineUpdateInterval';
end
data = load(file_path, load_names{:});
trigger_history = data.(field_name);
if isempty(trigger_history) || ndims(trigger_history) ~= 3
    return;
end

agent_count = size(trigger_history, 1);
iter_count = size(trigger_history, 2);
sample_count = size(trigger_history, 3);
time = data.t_set(:).';
time = time(1:sample_count);
if isfield(data, 'IPOnlineUpdateInterval') && data.IPOnlineUpdateInterval > 0
    inner_time_span = data.IPOnlineUpdateInterval;
else
    inner_time_span = median(diff(data.t_set));
end
inner_dt = inner_time_span / max(iter_count, 1);

event_times = cell(agent_count, 1);
event_counts = zeros(agent_count, 1);
for agent_i = 1:agent_count
    agent_history = squeeze(trigger_history(agent_i,:,:));
    [iter_idx, sample_idx] = find(agent_history);
    event_times{agent_i} = time(sample_idx(:)).' + (iter_idx(:).' - 1) * inner_dt;
    event_counts(agent_i) = numel(event_times{agent_i});
end

history.event_times = event_times;
history.event_counts = event_counts;
history.agent_count = agent_count;
history.time = data.t_set(:).';
history.source = file_path;
history.kind = 'inner';
end

function history = load_trigger_history_from_diagnostics( ...
    diagnostics_folder, method, form_tag, field_name, allow_ac_file, use_inner_history)
history = empty_trigger_history();
if ~isfolder(diagnostics_folder)
    return;
end

pattern = sprintf('diag_%s*%s*.mat', lower(method), form_tag);
candidates = [
    dir(fullfile(diagnostics_folder, pattern));
    dir(fullfile(diagnostics_folder, sprintf('diag_%s*.mat', lower(method))))];
if isempty(candidates)
    return;
end

[~, order] = sort([candidates.datenum], 'descend');
candidates = candidates(order);
for idx = 1:numel(candidates)
    if ~allow_ac_file && contains(candidates(idx).name, '_ac_')
        continue;
    end
    candidate_file = fullfile(candidates(idx).folder, candidates(idx).name);
    if use_inner_history
        history = load_inner_trigger_history(candidate_file, field_name);
    else
        history = load_trigger_history(candidate_file, field_name);
    end
    if ~isempty(history.counts) || ~isempty(history.event_times)
        fprintf('Using %s from diagnostics file: %s\n', field_name, candidate_file);
        return;
    end
end

history = empty_trigger_history();
end

function plot_trigger_instances(ax, history, label_prefix)
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
set(ax, 'FontName', 'Times New Roman', 'FontSize', 11);
if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
    ax.Toolbar.Visible = 'off';
end

if isempty(history.counts) && isempty(history.event_times)
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Agent');
    text(ax, 0.5, 0.5, sprintf('No %s data', label_prefix), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
    return;
end

if ~isempty(history.event_times)
    agent_count = history.agent_count;
else
    agent_count = size(history.counts, 1);
end
agent_colors = lines(agent_count);
for agent_i = 1:agent_count
    if ~isempty(history.event_times)
        trigger_x = history.event_times{agent_i};
        total_broadcast_count = history.event_counts(agent_i);
    else
        event_idx = find(history.counts(agent_i, :) > 0);
        trigger_x = history.time(event_idx);
        total_broadcast_count = sum(history.counts(agent_i, :));
    end
    if isempty(trigger_x)
        continue;
    end

    trigger_y = agent_i * ones(size(trigger_x));
    plot(ax, trigger_x, trigger_y, '*', ...
        'Color', agent_colors(agent_i,:), 'MarkerSize', 5.5, ...
        'DisplayName', sprintf('Agent %d, n=%.0f', agent_i, total_broadcast_count));
end

if ~isempty(history.event_times)
    total_event_count = sum(history.event_counts);
else
    total_event_count = nnz(history.counts > 0);
end
if total_event_count == 0
    text(ax, 0.5, 0.5, sprintf('No %s broadcast events recorded', label_prefix), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
end

xlabel(ax, 'Time (s)');
ylabel(ax, 'Agent');
ylim(ax, [0.5, agent_count + 0.5]);
yticks(ax, 1:agent_count);
if ~isempty(history.time)
    xlim(ax, [history.time(1), history.time(end)]);
end
if total_event_count > 0
    legend(ax, 'Location', 'eastoutside', 'FontSize', 8);
end
end

function history = empty_trigger_history()
history = struct('counts', [], 'time', [], 'source', '', ...
    'event_times', [], 'event_counts', [], 'agent_count', [], 'kind', '');
end
