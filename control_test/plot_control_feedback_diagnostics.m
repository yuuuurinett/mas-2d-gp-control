function plot_control_feedback_diagnostics(result_folder, output_folder)
%PLOT_CONTROL_FEEDBACK_DIAGNOSTICS Advisor-facing control diagnostics.
%
% Generates:
%   1) no-GP divergence plot for the baseline with unknown dynamics retained;
%   2) IP-DAC MC10 binned broadcast counts over physical simulation time.

if nargin < 1 || isempty(result_folder)
    result_folder = fullfile('control_test', 'result', 'Control_MC10');
end
if nargin < 2 || isempty(output_folder)
    output_folder = fullfile(result_folder, 'Figures');
end
if ~isfolder(output_folder)
    mkdir(output_folder);
end

diagnostics_folder = fullfile('control_test', 'result', 'Diagnostics');
plot_without_gp_divergence(diagnostics_folder, output_folder);
plot_dac_binned_triggers(result_folder, output_folder);
end

function plot_without_gp_divergence(diagnostics_folder, output_folder)
candidate = fullfile(diagnostics_folder, ...
    'diag_without_gp_unknown_on_formation_seed1001.mat');
if ~isfile(candidate)
    warning('No without-GP diagnostic file found: %s', candidate);
    return;
end

d = load(candidate, 'TrackingError_vector', 't_set', ...
    'vartheta_all_set', 'online_update_count');
t = d.t_set(:).';
agent_error = compute_agent_tracking_error(d);
sample_count = min(numel(t), size(agent_error, 2));
t = t(1:sample_count);
agent_error = agent_error(:,1:sample_count);

finite_mask = isfinite(agent_error);
first_nonfinite_idx = find(any(~finite_mask, 1), 1, 'first');

fig = figure('Name', 'Without GP divergence', 'Color', 'w', ...
    'Position', [100 100 820 470]);
ax = axes(fig);
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
set(ax, 'FontName', 'Times New Roman', 'FontSize', 12);

agent_count = size(agent_error, 1);
agent_colors = lines(agent_count);
for agent_i = 1:agent_count
    valid_idx = isfinite(agent_error(agent_i,:));
    plot(ax, t(valid_idx), agent_error(agent_i,valid_idx), '-*', ...
        'Color', agent_colors(agent_i,:), 'LineWidth', 1.2, ...
        'MarkerSize', 4.5, 'MarkerIndices', 1:20:nnz(valid_idx), ...
        'DisplayName', sprintf('Agent %d', agent_i));
end
if ~isempty(first_nonfinite_idx)
    xline(ax, t(first_nonfinite_idx), 'r--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('first nonfinite: %.2f s', t(first_nonfinite_idx)));
end

xlabel(ax, 'Time (s)');
ylabel(ax, 'Per-agent tracking error ||\vartheta_i(t)||');
title(ax, 'Without GP divergence under unknown dynamics', ...
    'FontWeight', 'bold');
legend(ax, 'Location', 'eastoutside');

saveas(fig, fullfile(output_folder, 'WithoutGP_unknown_on_agent_divergence.png'));
savefig(fig, fullfile(output_folder, 'WithoutGP_unknown_on_agent_divergence.fig'));
end

function agent_error = compute_agent_tracking_error(d)
if ~isfield(d, 'vartheta_all_set')
    agent_error = d.TrackingError_vector(:).';
    return;
end

if isfield(d, 'online_update_count') && ~isempty(d.online_update_count)
    agent_count = numel(d.online_update_count);
else
    agent_count = 6;
end

state_count = size(d.vartheta_all_set, 1);
if mod(state_count, agent_count) ~= 0
    agent_error = d.TrackingError_vector(:).';
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

function plot_dac_binned_triggers(result_folder, output_folder)
methods = {'poe','gpoe','moe','bcm','rbcm'};
labels = {'POE','GPOE','MOE','BCM','RBCM'};
bin_edges = 0:1:10;
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
bin_values = nan(numel(methods), numel(bin_centers));

for method_idx = 1:numel(methods)
    pattern = sprintf('%s_formation_mc*_seed*.mat', methods{method_idx});
    files = dir(fullfile(result_folder, pattern));
    if isempty(files)
        continue;
    end

    trial_counts = [];
    t_ref = [];
    for file_idx = 1:numel(files)
        d = load(fullfile(files(file_idx).folder, files(file_idx).name), ...
            'dac_broadcast_count_set', 't_set');
        if ~isfield(d, 'dac_broadcast_count_set') || ~isfield(d, 't_set')
            continue;
        end
        counts_per_step = mean(d.dac_broadcast_count_set, 1);
        trial_counts(file_idx, 1:numel(counts_per_step)) = counts_per_step; %#ok<AGROW>
        t_ref = d.t_set(1:numel(counts_per_step));
    end
    if isempty(trial_counts)
        continue;
    end

    mc_mean_counts = mean(trial_counts, 1, 'omitnan');
    for bin_idx = 1:numel(bin_centers)
        mask = t_ref >= bin_edges(bin_idx) & t_ref < bin_edges(bin_idx+1);
        bin_values(method_idx, bin_idx) = sum(mc_mean_counts(mask));
    end
end

fig = figure('Name', 'IP-DAC binned trigger counts MC10', 'Color', 'w', ...
    'Position', [100 100 840 460]);
ax = axes(fig);
bar(ax, bin_centers, bin_values.', 'grouped');
grid(ax, 'on'); box(ax, 'on');
set(ax, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel(ax, 'Time (s)');
ylabel(ax, 'Avg. broadcasts / agent per 1-s bin');
title(ax, 'IP-DAC average broadcast trigger times (Monte Carlo 10)', ...
    'FontWeight', 'bold');
legend(ax, labels, 'Location', 'northoutside', 'Orientation', 'horizontal');
xlim(ax, [bin_edges(1), bin_edges(end)]);

saveas(fig, fullfile(output_folder, 'IP_DAC_MC10_binned_broadcasts.png'));
savefig(fig, fullfile(output_folder, 'IP_DAC_MC10_binned_broadcasts.fig'));
end
