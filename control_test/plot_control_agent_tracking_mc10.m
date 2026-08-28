function plot_control_agent_tracking_mc10(result_folder, output_folder, method)
%PLOT_CONTROL_AGENT_TRACKING_MC10 Plot MC tracking-error curves per agent.
%
% The figure follows the advisor's shared-code convention:
%   e_i(t) = ||vartheta_i(t)||
% is computed for every agent and every Monte Carlo trial, then averaged
% across trials at each physical simulation time.

if nargin < 1 || isempty(result_folder)
    script_folder = fileparts(mfilename('fullpath'));
    project_root = fileparts(script_folder);
    result_folder = choose_default_result_folder(project_root, script_folder);
end
if nargin < 2 || isempty(output_folder)
    output_folder = fullfile(result_folder, 'Figures');
end
if nargin < 3 || isempty(method)
    method = 'poe';
end
if ~isfolder(output_folder)
    mkdir(output_folder);
end

method = lower(string(method));
mode_specs = [
    make_mode(method, "IP-DAC online", [0.0000 0.4470 0.7410])
    make_mode(method + "_ac", "IP-AC online", [0.8500 0.3250 0.0980])
    make_mode(method + "_offline", "IP-DAC offline", [0.4660 0.6740 0.1880])
    make_mode(method + "_ac_offline", "IP-AC offline", [0.4940 0.1840 0.5560])
    ];

[t_ref, agent_count] = find_reference_time_and_agent_count(result_folder, mode_specs);
if isempty(t_ref)
    warning('No MC10 control result files found in %s.', result_folder);
    return;
end

plot_system_average_convergence(result_folder, output_folder, method, ...
    mode_specs, t_ref, agent_count);

fig = figure('Name', sprintf('%s MC10 per-agent tracking error', upper(method)), ...
    'Color', 'w', 'Position', [80 60 1050 760]);
layout = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for agent_idx = 1:agent_count
    ax = nexttile(layout);
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 11, ...
        'YScale', 'log');

    for mode_idx = 1:numel(mode_specs)
        [mean_curve, min_curve, max_curve] = load_mode_agent_curve( ...
            result_folder, mode_specs(mode_idx).key, t_ref, agent_idx);
        if isempty(mean_curve)
            continue;
        end
        color_i = mode_specs(mode_idx).color;
        fill(ax, [t_ref, fliplr(t_ref)], [min_curve, fliplr(max_curve)], ...
            color_i, 'FaceAlpha', 0.16, 'EdgeAlpha', 0, ...
            'HandleVisibility', 'off');
        plot(ax, t_ref, mean_curve, '-', 'Color', color_i, ...
            'LineWidth', 1.4, 'DisplayName', mode_specs(mode_idx).label);
    end

    title(ax, sprintf('Agent %d', agent_idx));
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'e_i(t)');
    ylim(ax, [1e-3, 2]);
    if agent_idx == 1
        legend(ax, 'Location', 'best');
    end
end

title(layout, sprintf('%s per-agent tracking error, Monte Carlo 10', ...
    upper(method)), 'FontName', 'Times New Roman', ...
    'FontWeight', 'bold');

png_file = fullfile(output_folder, sprintf( ...
    'Control_MC10_%s_per_agent_tracking_error.png', method));
fig_file = fullfile(output_folder, sprintf( ...
    'Control_MC10_%s_per_agent_tracking_error.fig', method));
saveas(fig, png_file);
savefig(fig, fig_file);
fprintf('Saved per-agent tracking-error figure:\n%s\n%s\n', png_file, fig_file);
end

function result_folder = choose_default_result_folder(project_root, script_folder)
candidate_names = {'Control_MC10_T10_M800_th035_rbcm075', ...
    'Control_MC10_T10_M800_th035', 'Control_MC10_T20', 'Control_MC10'};
for name_idx = 1:numel(candidate_names)
    folder_name = candidate_names{name_idx};
    project_result_folder = fullfile(project_root, 'result', folder_name);
    legacy_result_folder = fullfile(script_folder, 'result', folder_name);
    if isfolder(legacy_result_folder) && ...
            ~isempty(dir(fullfile(legacy_result_folder, '*.mat')))
        result_folder = legacy_result_folder;
        return;
    end
    if isfolder(project_result_folder) && ...
            ~isempty(dir(fullfile(project_result_folder, '*.mat')))
        result_folder = project_result_folder;
        return;
    end
end
result_folder = fullfile(script_folder, 'result', candidate_names{1});
end

function plot_system_average_convergence(result_folder, output_folder, method, ...
    mode_specs, t_ref, agent_count)
fig = figure('Name', sprintf('%s MC10 average tracking error', upper(method)), ...
    'Color', 'w', 'Position', [120 120 720 480]);
ax = axes(fig);
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
set(ax, 'FontName', 'Times New Roman', 'FontSize', 12, 'YScale', 'log');

for mode_idx = 1:numel(mode_specs)
    [mean_curve, min_curve, max_curve] = load_mode_average_curve( ...
        result_folder, mode_specs(mode_idx).key, t_ref, agent_count);
    if isempty(mean_curve)
        continue;
    end
    color_i = mode_specs(mode_idx).color;
    fill(ax, [t_ref, fliplr(t_ref)], [min_curve, fliplr(max_curve)], ...
        color_i, 'FaceAlpha', 0.16, 'EdgeAlpha', 0, ...
        'HandleVisibility', 'off');
    plot(ax, t_ref, mean_curve, '-', 'Color', color_i, ...
        'LineWidth', 1.8, 'DisplayName', mode_specs(mode_idx).label);
end

title(ax, sprintf('%s average tracking error, Monte Carlo 10', upper(method)), ...
    'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel(ax, 'Time (s)');
ylabel(ax, 'Agent-average tracking error');
ylim(ax, [1e-3, 2]);
legend(ax, 'Location', 'best');

png_file = fullfile(output_folder, sprintf( ...
    'Control_MC10_%s_average_tracking_error.png', method));
fig_file = fullfile(output_folder, sprintf( ...
    'Control_MC10_%s_average_tracking_error.fig', method));
saveas(fig, png_file);
savefig(fig, fig_file);
fprintf('Saved average tracking-error figure:\n%s\n%s\n', png_file, fig_file);
end

function [mean_curve, min_curve, max_curve] = load_mode_average_curve( ...
    result_folder, mode_key, t_ref, agent_count)
files = dir(fullfile(result_folder, sprintf('%s_formation_mc*_seed*.mat', ...
    mode_key)));
if isempty(files)
    mean_curve = [];
    min_curve = [];
    max_curve = [];
    return;
end

trial_curve_set = nan(numel(files), numel(t_ref));
for file_idx = 1:numel(files)
    d = load(fullfile(files(file_idx).folder, files(file_idx).name), ...
        't_set', 'vartheta_all_set', 'online_update_count');
    agent_error = compute_agent_tracking_error(d);
    if isempty(agent_error)
        continue;
    end
    agent_count_i = min(agent_count, size(agent_error, 1));
    mean_agent_error = mean(agent_error(1:agent_count_i,:), 1, 'omitnan');
    t = d.t_set(:).';
    finite_mask = isfinite(mean_agent_error) & isfinite(t);
    if nnz(finite_mask) < 2
        continue;
    end
    trial_curve_set(file_idx,:) = interp1(t(finite_mask), ...
        mean_agent_error(finite_mask), t_ref, 'linear', NaN);
end

if all(isnan(trial_curve_set), 'all')
    mean_curve = [];
    min_curve = [];
    max_curve = [];
    return;
end

mean_curve = mean(trial_curve_set, 1, 'omitnan');
std_curve = std(trial_curve_set, 0, 1, 'omitnan');
min_curve = max(mean_curve - std_curve, 1e-10);
max_curve = mean_curve + std_curve;
end

function mode_spec = make_mode(key, label, color)
mode_spec = struct('key', string(key), 'label', string(label), ...
    'color', color);
end

function [t_ref, agent_count] = find_reference_time_and_agent_count( ...
    result_folder, mode_specs)
t_ref = [];
agent_count = [];
for mode_idx = 1:numel(mode_specs)
    files = dir(fullfile(result_folder, sprintf('%s_formation_mc*_seed*.mat', ...
        mode_specs(mode_idx).key)));
    for file_idx = 1:numel(files)
        d = load(fullfile(files(file_idx).folder, files(file_idx).name), ...
            't_set', 'vartheta_all_set', 'online_update_count');
        if isfield(d, 't_set') && isfield(d, 'vartheta_all_set')
            t_ref = d.t_set(:).';
            agent_count = infer_agent_count(d);
            return;
        end
    end
end
end

function [mean_curve, min_curve, max_curve] = load_mode_agent_curve( ...
    result_folder, mode_key, t_ref, agent_idx)
files = dir(fullfile(result_folder, sprintf('%s_formation_mc*_seed*.mat', ...
    mode_key)));
if isempty(files)
    mean_curve = [];
    min_curve = [];
    max_curve = [];
    return;
end

trial_curve_set = nan(numel(files), numel(t_ref));
for file_idx = 1:numel(files)
    d = load(fullfile(files(file_idx).folder, files(file_idx).name), ...
        't_set', 'vartheta_all_set', 'online_update_count');
    agent_error = compute_agent_tracking_error(d);
    if isempty(agent_error) || agent_idx > size(agent_error, 1)
        continue;
    end
    t = d.t_set(:).';
    finite_mask = isfinite(agent_error(agent_idx,:)) & isfinite(t);
    if nnz(finite_mask) < 2
        continue;
    end
    trial_curve_set(file_idx,:) = interp1(t(finite_mask), ...
        agent_error(agent_idx,finite_mask), t_ref, 'linear', NaN);
end

if all(isnan(trial_curve_set), 'all')
    mean_curve = [];
    min_curve = [];
    max_curve = [];
    return;
end

mean_curve = mean(trial_curve_set, 1, 'omitnan');
std_curve = std(trial_curve_set, 0, 1, 'omitnan');
min_curve = max(mean_curve - std_curve, 1e-10);
max_curve = mean_curve + std_curve;
end

function agent_error = compute_agent_tracking_error(d)
agent_error = [];
if ~isfield(d, 'vartheta_all_set')
    return;
end

agent_count = infer_agent_count(d);
state_count = size(d.vartheta_all_set, 1);
if agent_count <= 0 || mod(state_count, agent_count) ~= 0
    return;
end

state_dim = state_count / agent_count;
sample_count = size(d.vartheta_all_set, 2);
agent_error = nan(agent_count, sample_count);
for sample_idx = 1:sample_count
    vartheta_sample = reshape(d.vartheta_all_set(:,sample_idx), ...
        state_dim, agent_count);
    agent_error(:,sample_idx) = vecnorm(vartheta_sample, 2, 1).';
end
end

function agent_count = infer_agent_count(d)
if isfield(d, 'online_update_count') && ~isempty(d.online_update_count)
    agent_count = numel(d.online_update_count);
else
    agent_count = 6;
end
end
