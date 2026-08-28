function saved_files = plot_prediction_true_dynamic_mc10(result_folder, output_folder, method, mode_key, mc_index)
%PLOT_PREDICTION_TRUE_DYNAMIC_MC10 Plot GP prediction vs true dynamics for MC control runs.
%
% For each agent, this function loads all Monte Carlo trials for one control
% mode, averages f_true_all_set and f_hat_all_set across trials at each
% physical simulation time, and saves one figure per agent.
%
% Default mode is online IP-DAC for the selected aggregation method. If
% mc_index is provided, one Monte Carlo trial is plotted directly instead
% of averaging trajectories from different initial conditions.

if nargin < 1 || isempty(result_folder)
    script_folder = fileparts(mfilename('fullpath'));
    project_root = fileparts(script_folder);
    result_folder = choose_default_result_folder(project_root);
end
if nargin < 2 || isempty(output_folder)
    output_folder = fullfile(result_folder, 'Figures', 'PredictionTrue');
end
if nargin < 3 || isempty(method)
    method = 'poe';
end
if nargin < 4 || isempty(mode_key)
    mode_key = lower(string(method));
else
    mode_key = lower(string(mode_key));
end
if nargin < 5
    mc_index = [];
end

method = lower(string(method));
mode_label = make_mode_label(method, mode_key);
if ~isfolder(output_folder)
    mkdir(output_folder);
end

files = dir(fullfile(result_folder, sprintf('%s_formation_mc*_seed*.mat', mode_key)));
if isempty(files)
    warning('No files found for mode "%s" in %s.', mode_key, result_folder);
    saved_files = strings(0,1);
    return;
end
files = filter_mc_files(files, mc_index);
if isempty(files)
    warning('No files found for mode "%s" and requested MC index.', mode_key);
    saved_files = strings(0,1);
    return;
end

[t_ref, y_dim, agent_count] = find_reference_metadata(files);
if isempty(t_ref)
    warning('No valid f_true_all_set/f_hat_all_set data found for mode "%s".', mode_key);
    saved_files = strings(0,1);
    return;
end

[f_true_mean, f_true_std, f_hat_mean, f_hat_std] = load_mc_prediction_curves( ...
    files, t_ref, y_dim, agent_count);
is_single_trial = isscalar(files);
if is_single_trial
    trial_label = extract_trial_label(files(1).name);
else
    trial_label = sprintf('Monte Carlo %d', numel(files));
end

saved_files = strings(agent_count*2, 1);
save_idx = 0;
for agent_idx = 1:agent_count
    fig = figure('Name', sprintf('%s agent %d prediction vs true', upper(mode_key), agent_idx), ...
        'Color', 'w', 'Position', [100 80 860 560]);
    layout = tiledlayout(fig, y_dim, 1, 'TileSpacing', 'compact', ...
        'Padding', 'compact');

    for dim_idx = 1:y_dim
        ax = nexttile(layout);
        hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
        set(ax, 'FontName', 'Times New Roman', 'FontSize', 11);

        true_curve = squeeze(f_true_mean(dim_idx, agent_idx, :)).';
        pred_curve = squeeze(f_hat_mean(dim_idx, agent_idx, :)).';
        true_std_curve = squeeze(f_true_std(dim_idx, agent_idx, :)).';
        pred_std_curve = squeeze(f_hat_std(dim_idx, agent_idx, :)).';

        if is_single_trial
            plot(ax, t_ref, true_curve, '-', 'Color', [0.15 0.15 0.15], ...
                'LineWidth', 1.8, 'DisplayName', 'True dynamic');
            plot(ax, t_ref, pred_curve, '-', 'Color', [0.0000 0.4470 0.7410], ...
                'LineWidth', 1.5, 'DisplayName', 'GP prediction');
        else
            plot_shaded_curve(ax, t_ref, true_curve, true_std_curve, ...
                [0.15 0.15 0.15], 'True dynamic');
            plot_shaded_curve(ax, t_ref, pred_curve, pred_std_curve, ...
                [0.0000 0.4470 0.7410], 'GP prediction');
        end

        ylabel(ax, sprintf('f_{%d,%d}(t)', agent_idx, dim_idx));
        if dim_idx == y_dim
            xlabel(ax, 'Time (s)');
        end
        if dim_idx == 1
            legend(ax, 'Location', 'best');
        end
        xlim(ax, [t_ref(1), t_ref(end)]);
    end

    title(layout, sprintf('%s %s: prediction vs true dynamic, %s, Agent %d', ...
        upper(method), mode_label, trial_label, agent_idx), ...
        'FontName', 'Times New Roman', 'FontWeight', 'bold');

    if is_single_trial
        file_tag = sprintf('Control_%s_%s_prediction_true_agent%d', ...
            trial_label, mode_key, agent_idx);
        file_tag = regexprep(file_tag, '[^A-Za-z0-9_]', '_');
    else
        file_tag = sprintf('Control_MC%d_%s_prediction_true_agent%d', ...
            numel(files), mode_key, agent_idx);
    end
    png_file = fullfile(output_folder, sprintf( ...
        '%s.png', file_tag));
    fig_file = fullfile(output_folder, sprintf( ...
        '%s.fig', file_tag));
    saveas(fig, png_file);
    savefig(fig, fig_file);
    close(fig);

    save_idx = save_idx + 1;
    saved_files(save_idx) = string(png_file);
    save_idx = save_idx + 1;
    saved_files(save_idx) = string(fig_file);
end

saved_files = saved_files(1:save_idx);
fprintf('Saved %d prediction-vs-true files for %s in:\n%s\n', ...
    numel(saved_files), mode_key, output_folder);
end

function files = filter_mc_files(files, mc_index)
if isempty(mc_index)
    return;
end
mc_index = double(mc_index);
keep_mask = false(numel(files), 1);
for file_idx = 1:numel(files)
    token = regexp(files(file_idx).name, '_mc(\d+)_seed', 'tokens', 'once');
    keep_mask(file_idx) = ~isempty(token) && str2double(token{1}) == mc_index;
end
files = files(keep_mask);
end

function trial_label = extract_trial_label(file_name)
token = regexp(file_name, '_mc(\d+)_seed(\d+)', 'tokens', 'once');
if isempty(token)
    trial_label = 'single trial';
else
    trial_label = sprintf('MC%s seed%s', token{1}, token{2});
end
end

function mode_label = make_mode_label(method, mode_key)
method = lower(string(method));
mode_key = lower(string(mode_key));
if mode_key == method
    mode_label = 'IP-DAC online';
elseif mode_key == method + "_ac"
    mode_label = 'IP-AC online';
elseif mode_key == method + "_offline"
    mode_label = 'IP-DAC offline';
elseif mode_key == method + "_ac_offline"
    mode_label = 'IP-AC offline';
else
    mode_label = char(upper(mode_key));
end
end

function result_folder = choose_default_result_folder(project_root)
candidate_names = {'Control_MC10_T10_M800_th035_rbcm075', ...
    'Control_MC10_T10_M800_th035', 'Control_MC10_T20', 'Control_MC10'};
result_folder = fullfile(project_root, 'result', candidate_names{1});
for name_idx = 1:numel(candidate_names)
    candidate_folder = fullfile(project_root, 'result', candidate_names{name_idx});
    if isfolder(candidate_folder) && ~isempty(dir(fullfile(candidate_folder, '*.mat')))
        result_folder = candidate_folder;
        return;
    end
end
end

function [t_ref, y_dim, agent_count] = find_reference_metadata(files)
t_ref = [];
y_dim = [];
agent_count = [];
for file_idx = 1:numel(files)
    d = load(fullfile(files(file_idx).folder, files(file_idx).name), ...
        't_set', 'f_true_all_set', 'f_hat_all_set');
    if isfield(d, 't_set') && isfield(d, 'f_true_all_set') && ...
            isfield(d, 'f_hat_all_set')
        t_ref = d.t_set(:).';
        y_dim = size(d.f_true_all_set, 1);
        agent_count = size(d.f_true_all_set, 2);
        return;
    end
end
end

function [f_true_mean, f_true_std, f_hat_mean, f_hat_std] = load_mc_prediction_curves( ...
    files, t_ref, y_dim, agent_count)
trial_count = numel(files);
sample_count = numel(t_ref);
f_true_trials = nan(y_dim, agent_count, sample_count, trial_count);
f_hat_trials = nan(y_dim, agent_count, sample_count, trial_count);

for file_idx = 1:trial_count
    d = load(fullfile(files(file_idx).folder, files(file_idx).name), ...
        't_set', 'f_true_all_set', 'f_hat_all_set');
    if ~isfield(d, 't_set') || ~isfield(d, 'f_true_all_set') || ...
            ~isfield(d, 'f_hat_all_set')
        continue;
    end
    [true_interp, hat_interp] = interpolate_prediction_set(d, t_ref, y_dim, agent_count);
    f_true_trials(:,:,:,file_idx) = true_interp;
    f_hat_trials(:,:,:,file_idx) = hat_interp;
end

f_true_mean = mean(f_true_trials, 4, 'omitnan');
f_true_std = std(f_true_trials, 0, 4, 'omitnan');
f_hat_mean = mean(f_hat_trials, 4, 'omitnan');
f_hat_std = std(f_hat_trials, 0, 4, 'omitnan');
end

function [true_interp, hat_interp] = interpolate_prediction_set(d, t_ref, y_dim, agent_count)
sample_count = numel(t_ref);
true_interp = nan(y_dim, agent_count, sample_count);
hat_interp = nan(y_dim, agent_count, sample_count);
t = d.t_set(:).';
for dim_idx = 1:y_dim
    for agent_idx = 1:agent_count
        true_curve = squeeze(d.f_true_all_set(dim_idx, agent_idx, :)).';
        hat_curve = squeeze(d.f_hat_all_set(dim_idx, agent_idx, :)).';
        true_interp(dim_idx, agent_idx, :) = interpolate_curve(t, true_curve, t_ref);
        hat_interp(dim_idx, agent_idx, :) = interpolate_curve(t, hat_curve, t_ref);
    end
end
end

function y_ref = interpolate_curve(t, y, t_ref)
finite_mask = isfinite(t) & isfinite(y);
if nnz(finite_mask) < 2
    y_ref = nan(1, numel(t_ref));
    return;
end
y_ref = interp1(t(finite_mask), y(finite_mask), t_ref, 'linear', NaN);
end

function plot_shaded_curve(ax, t, mean_curve, std_curve, color_i, label_text)
lower_curve = mean_curve - std_curve;
upper_curve = mean_curve + std_curve;
fill(ax, [t, fliplr(t)], [lower_curve, fliplr(upper_curve)], ...
    color_i, 'FaceAlpha', 0.14, 'EdgeAlpha', 0, ...
    'HandleVisibility', 'off');
plot(ax, t, mean_curve, '-', 'Color', color_i, ...
    'LineWidth', 1.6, 'DisplayName', label_text);
end
