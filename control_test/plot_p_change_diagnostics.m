function plot_p_change_diagnostics(result_file, output_folder)
%PLOT_P_CHANGE_DIAGNOSTICS Plot GP inducing-information changes vs trigger threshold.
%
% The physical-time broadcast trigger in run_simulation_inducing_point uses
% the RMS change of each agent's GP inducing-point information P relative to
% its last broadcast reference.  This diagnostic separates total P-change
% into numerator-like and denominator/precision-like components.

if nargin < 1 || isempty(result_file)
    script_folder = fileparts(mfilename('fullpath'));
    result_folder = fullfile(script_folder, 'result', 'Diagnostics');
    candidates = dir(fullfile(result_folder, 'diag_poe_M800_Pchange_T10_seed*.mat'));
    if isempty(candidates)
        candidates = dir(fullfile(result_folder, 'diag_poe*Pchange*.mat'));
    end
    if isempty(candidates)
        error('No P-change diagnostic .mat file found in %s.', result_folder);
    end
    [~, order] = sort([candidates.datenum], 'descend');
    result_file = fullfile(candidates(order(1)).folder, candidates(order(1)).name);
end

if nargin < 2 || isempty(output_folder)
    script_folder = fileparts(mfilename('fullpath'));
    output_folder = fullfile(script_folder, 'result', 'Figures');
end
if ~isfolder(output_folder)
    mkdir(output_folder);
end

required_fields = ["t_set", ...
    "consensus_input_change_per_agent_set", ...
    "consensus_input_numerator_change_per_agent_set", ...
    "consensus_input_denominator_change_per_agent_set", ...
    "ConsensusInputTriggerTol"];
names = string({whos('-file', result_file).name});
missing_fields = required_fields(~ismember(required_fields, names));
if ~isempty(missing_fields)
    error('Missing required fields in %s: %s', result_file, strjoin(missing_fields, ', '));
end

data = load(result_file, required_fields{:}, ...
    'dac_broadcast_count_set', 'ac_physical_trigger_count_set', ...
    'CurrentMode', 'use_formation');

t = data.t_set(:).';
sample_count = size(data.consensus_input_change_per_agent_set, 2);
t = t(1:sample_count);
threshold = data.ConsensusInputTriggerTol;

total_change = data.consensus_input_change_per_agent_set;
num_change = data.consensus_input_numerator_change_per_agent_set;
den_change = data.consensus_input_denominator_change_per_agent_set;

finite_total = total_change;
finite_num = num_change;
finite_den = den_change;
finite_total(~isfinite(finite_total)) = NaN;
finite_num(~isfinite(finite_num)) = NaN;
finite_den(~isfinite(finite_den)) = NaN;

valid_time_idx = any(isfinite(finite_total), 1);
t_change = t(valid_time_idx);
finite_total = finite_total(:,valid_time_idx);
finite_num = finite_num(:,valid_time_idx);
finite_den = finite_den(:,valid_time_idx);

total_median = median(finite_total, 1, 'omitnan');
num_median = median(finite_num, 1, 'omitnan');
den_median = median(finite_den, 1, 'omitnan');
total_max = max(finite_total, [], 1, 'omitnan');

trigger_counts = [];
trigger_label = 'trigger';
if isfield(data, 'dac_broadcast_count_set') && ~isempty(data.dac_broadcast_count_set)
    trigger_counts = data.dac_broadcast_count_set;
    trigger_label = 'DAC broadcast';
elseif isfield(data, 'ac_physical_trigger_count_set') && ~isempty(data.ac_physical_trigger_count_set)
    trigger_counts = data.ac_physical_trigger_count_set;
    trigger_label = 'AC broadcast';
end

fig = figure('Name', 'P-change diagnostics', 'Color', 'w', ...
    'Position', [100 100 980 620]);
layout = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(layout);
hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');
set(ax1, 'FontName', 'Times New Roman', 'FontSize', 11);
if isprop(ax1, 'Toolbar') && ~isempty(ax1.Toolbar)
    ax1.Toolbar.Visible = 'off';
end
plot(ax1, t_change, total_median, 'k-o', 'LineWidth', 1.5, ...
    'MarkerSize', 3.0, ...
    'DisplayName', 'median total P-change');
plot(ax1, t_change, num_median, 'b-o', 'LineWidth', 1.2, ...
    'MarkerSize', 3.0, ...
    'DisplayName', 'median numerator change');
plot(ax1, t_change, den_median, 'r-o', 'LineWidth', 1.2, ...
    'MarkerSize', 3.0, ...
    'DisplayName', 'median denominator change');
plot(ax1, [t(1), t(end)], threshold * [1, 1], 'k--', ...
    'LineWidth', 1.0, 'DisplayName', sprintf('threshold = %.2f', threshold));
ylabel(ax1, 'RMS change');
title(ax1, 'Median P-change across agents', 'Interpreter', 'none');
legend(ax1, 'Location', 'eastoutside', 'FontSize', 8);

ax2 = nexttile(layout);
hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
set(ax2, 'FontName', 'Times New Roman', 'FontSize', 11);
if isprop(ax2, 'Toolbar') && ~isempty(ax2.Toolbar)
    ax2.Toolbar.Visible = 'off';
end
plot(ax2, t_change, total_max, 'Color', [0.2 0.2 0.2], 'LineWidth', 1.4, ...
    'DisplayName', 'max total P-change');
plot(ax2, [t(1), t(end)], threshold * [1, 1], 'k--', ...
    'LineWidth', 1.0, 'DisplayName', sprintf('threshold = %.2f', threshold));
if ~isempty(trigger_counts)
    trigger_times = t(sum(trigger_counts, 1) > 0);
    trigger_y = threshold * ones(size(trigger_times));
    plot(ax2, trigger_times, trigger_y, '*', 'Color', [0 0.45 0.74], ...
        'MarkerSize', 4.5, 'DisplayName', trigger_label);
end
xlabel(ax2, 'Time (s)');
ylabel(ax2, 'RMS change');
title(ax2, 'Max P-change and physical-time broadcast instances', 'Interpreter', 'none');
legend(ax2, 'Location', 'eastoutside', 'FontSize', 8);

mode_label = 'unknown mode';
if isfield(data, 'CurrentMode')
    mode_label = string(data.CurrentMode);
end
form_label = 'formation';
if isfield(data, 'use_formation') && ~data.use_formation
    form_label = 'leader-only';
end
title(layout, sprintf('P-change diagnostics: %s [%s]', upper(mode_label), form_label), ...
    'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold', ...
    'Interpreter', 'none');

[~, base_name] = fileparts(result_file);
png_file = fullfile(output_folder, [base_name, '_Pchange.png']);
fig_file = fullfile(output_folder, [base_name, '_Pchange.fig']);
if isfile(png_file)
    delete(png_file);
end
exportgraphics(fig, png_file, 'Resolution', 200);
savefig(fig, fig_file);
fprintf('Saved P-change diagnostic figure: %s\n', png_file);
end
