function saved_files = plot_prediction_broadcast_diagnostics(result_file, output_folder)
%PLOT_PREDICTION_BROADCAST_DIAGNOSTICS Align prediction error with broadcasts.
%
% This single-run diagnostic is intended for advisor-facing debugging of
% event-triggered GP control. It avoids averaging f_true/f_hat trajectories
% from different Monte Carlo initial conditions.

if nargin < 1 || isempty(result_file)
    result_file = fullfile('result', 'Control_MC10_T10_M800_th035_rbcm075', ...
        'rbcm_formation_mc01_seed1001.mat');
end
if nargin < 2 || isempty(output_folder)
    [result_folder, result_name] = fileparts(result_file);
    output_folder = fullfile(result_folder, 'Figures', ...
        ['PredictionBroadcastDiagnostics_' result_name]);
end
if ~isfolder(output_folder)
    mkdir(output_folder);
end

d = load(result_file, 't_set', 'f_hat_all_set', 'f_true_all_set', ...
    'f_direct_all_set', 'dac_broadcast_count_set', 'ac_broadcast_count_set');
validate_required_fields(d, result_file);

t = d.t_set(:).';
f_hat = d.f_hat_all_set;
f_true = d.f_true_all_set;

sample_count = min([numel(t), size(f_hat, 3), size(f_true, 3)]);
t = t(1:sample_count);
f_hat = f_hat(:,:,1:sample_count);
f_true = f_true(:,:,1:sample_count);

prediction_error = squeeze(vecnorm(f_hat - f_true, 2, 1));
if isvector(prediction_error)
    prediction_error = prediction_error(:).';
end
has_direct_prediction = isfield(d, 'f_direct_all_set') && ...
    any(isfinite(d.f_direct_all_set(:)));
if has_direct_prediction
    f_direct = d.f_direct_all_set(:,:,1:sample_count);
    direct_prediction_error = squeeze(vecnorm(f_direct - f_true, 2, 1));
    projection_gap = squeeze(vecnorm(f_hat - f_direct, 2, 1));
    if isvector(direct_prediction_error)
        direct_prediction_error = direct_prediction_error(:).';
    end
    if isvector(projection_gap)
        projection_gap = projection_gap(:).';
    end
else
    direct_prediction_error = nan(size(prediction_error));
    projection_gap = nan(size(prediction_error));
end
prediction_jump = squeeze(vecnorm(diff(f_hat, 1, 3), 2, 1));
if isvector(prediction_jump)
    prediction_jump = prediction_jump(:).';
end

[broadcast_count_set, consensus_label] = choose_broadcast_count_set(d);
if isempty(broadcast_count_set)
    broadcast_count_set = zeros(size(prediction_error, 1), sample_count-1);
    consensus_label = 'no broadcast data';
end
broadcast_count_set = broadcast_count_set(:,1:min(size(broadcast_count_set,2), sample_count-1));

agent_count = size(prediction_error, 1);
fig = figure('Name', 'Prediction error and broadcast diagnostics', ...
    'Color', 'w', 'Position', [70 45 1080 760]);
layout = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

summary = table('Size', [agent_count, 5], ...
    'VariableTypes', {'double','double','double','double','double'}, ...
    'VariableNames', {'Agent','MeanPredictionError','MaxPredictionError', ...
    'MaxPredictionJump','BroadcastCount'});
if has_direct_prediction
    summary.MeanDirectPredictionError = nan(agent_count,1);
    summary.MeanProjectionGap = nan(agent_count,1);
end

for agent_idx = 1:agent_count
    ax = nexttile(layout);
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 11);

    err_i = prediction_error(agent_idx,:);
    direct_err_i = direct_prediction_error(agent_idx,:);
    gap_i = projection_gap(agent_idx,:);
    jump_i = prediction_jump(agent_idx,:);
    broadcast_i = broadcast_count_set(agent_idx,:) > 0;
    trigger_t = t(1:numel(broadcast_i));
    trigger_t = trigger_t(broadcast_i);

    plot(ax, t, err_i, '-', 'Color', [0.0000 0.4470 0.7410], ...
        'LineWidth', 1.5, 'DisplayName', 'masked ||f_{hat}-f_{true}||');
    if has_direct_prediction
        plot(ax, t, direct_err_i, '-', 'Color', [0.4660 0.6740 0.1880], ...
            'LineWidth', 1.25, 'DisplayName', 'direct ||f_{direct}-f_{true}||');
        plot(ax, t, gap_i, '--', 'Color', [0.4940 0.1840 0.5560], ...
            'LineWidth', 1.15, 'DisplayName', '||f_{hat}-f_{direct}||');
    end

    y_top = max([err_i(:); direct_err_i(:); gap_i(:)], [], 'omitnan');
    if ~isfinite(y_top) || y_top <= 0
        y_top = 1;
    end
    if ~isempty(trigger_t)
        plot(ax, trigger_t, y_top*ones(size(trigger_t)), 'r*', ...
            'MarkerSize', 5.5, 'LineWidth', 0.8, ...
            'DisplayName', 'broadcast');
    end

    [max_jump, jump_idx] = max(jump_i, [], 'omitnan');
    if ~isempty(jump_idx) && isfinite(max_jump) && jump_idx <= numel(t)
        plot(ax, t(jump_idx+1), err_i(jump_idx+1), 'ro', ...
            'MarkerSize', 6.5, 'LineWidth', 1.2, ...
            'DisplayName', 'max f_{hat} jump');
        text(ax, t(jump_idx+1), err_i(jump_idx+1), ...
            sprintf('  jump %.2f', max_jump), ...
            'FontName', 'Times New Roman', 'FontSize', 9, ...
            'Color', [0.65 0 0]);
    end

    title(ax, sprintf('Agent %d', agent_idx));
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Prediction error norm');
    xlim(ax, [t(1), t(end)]);
    ylim(ax, [0, y_top*1.18]);
    if agent_idx == 1
        legend(ax, 'Location', 'best');
    end

    summary.Agent(agent_idx) = agent_idx;
    summary.MeanPredictionError(agent_idx) = mean(err_i, 'omitnan');
    summary.MaxPredictionError(agent_idx) = max(err_i, [], 'omitnan');
    summary.MaxPredictionJump(agent_idx) = max_jump;
    summary.BroadcastCount(agent_idx) = sum(broadcast_i);
    if has_direct_prediction
        summary.MeanDirectPredictionError(agent_idx) = mean(direct_err_i, 'omitnan');
        summary.MeanProjectionGap(agent_idx) = mean(gap_i, 'omitnan');
    end
end

[~, result_name] = fileparts(result_file);
title(layout, sprintf('%s: prediction error and %s broadcasts', ...
    result_name, consensus_label), ...
    'Interpreter', 'none', 'FontName', 'Times New Roman', ...
    'FontWeight', 'bold');

png_file = fullfile(output_folder, [result_name '_prediction_error_broadcasts.png']);
fig_file = fullfile(output_folder, [result_name '_prediction_error_broadcasts.fig']);
csv_file = fullfile(output_folder, [result_name '_prediction_error_broadcasts_summary.csv']);
saveas(fig, png_file);
savefig(fig, fig_file);
writetable(summary, csv_file);

saved_files = string({png_file; fig_file; csv_file});
fprintf('Saved prediction-broadcast diagnostics:\n%s\n%s\n%s\n', ...
    png_file, fig_file, csv_file);
disp(summary);
if ~has_direct_prediction
    fprintf(['Note: f_direct_all_set is missing or all NaN. ', ...
        'Re-run the simulation with CONTROL_IP_PROJECTION_DIAGNOSTICS=1 ', ...
        'to separate direct local aggregation from masked-GP projection.\n']);
end
end

function validate_required_fields(d, result_file)
if ~isfield(d, 't_set') || ~isfield(d, 'f_hat_all_set') || ...
        ~isfield(d, 'f_true_all_set')
    error('Missing t_set/f_hat_all_set/f_true_all_set in %s.', result_file);
end
end

function [broadcast_count_set, consensus_label] = choose_broadcast_count_set(d)
if isfield(d, 'dac_broadcast_count_set') && ~isempty(d.dac_broadcast_count_set)
    broadcast_count_set = d.dac_broadcast_count_set;
    consensus_label = 'DAC';
elseif isfield(d, 'ac_broadcast_count_set') && ~isempty(d.ac_broadcast_count_set)
    broadcast_count_set = d.ac_broadcast_count_set;
    consensus_label = 'AC';
else
    broadcast_count_set = [];
    consensus_label = '';
end
end
