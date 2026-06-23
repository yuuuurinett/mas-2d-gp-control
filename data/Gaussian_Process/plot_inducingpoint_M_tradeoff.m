%% plot_inducingpoint_M_tradeoff_paperstyle.m
% Paper-style M-ablation summary:
%   Fixed setting: IP-DAC / event-triggered consensus
%   M = 100:100:2500
%   Metrics: SMSE and MSLL
%   Mean/std over Monte Carlo seeds
%
% Output:
%   1) Overall figure: M = 100:2500
%   2) Zoomed figure: M = 500:2500
%   3) CSV summary table

clear; clc; close all;

%% ===================== User settings =====================

Method_list  = {'poe','gpoe','moe','bcm','rbcm'};
Method_label = {'POE','GPOE','MOE','BCM','RBCM'};

Dataset_list = {'KIN40K','POL','PUMADYN32NM','SARCOS'};

M_list      = 100:100:2500;
seeds       = 1:3;
train_ratio = 0.4;
tr_tag      = round(train_ratio * 100);

% IMPORTANT:
% If this script is inside:
%   mas_2D_test/data/Gaussian_Process/
% and your full result folder is:
%   mas_2D_test/data/Gaussian_Process/Result/Dataset/...
% then this is correct.
ProjectRoot = fileparts(mfilename('fullpath'));

% If you want to force the path manually, uncomment this:
% ProjectRoot = '/Users/duyurou/Desktop/ip/mas_gp_code/dac_code/mas_2D_test/data/Gaussian_Process';

ResultRoot = fullfile(ProjectRoot, 'Result', 'Dataset');

FigFolder = fullfile(ProjectRoot, 'Result', 'Figures', 'M_ablation_paperstyle');
if ~exist(FigFolder, 'dir')
    mkdir(FigFolder);
end

show_errorbar = false;   % false = cleaner paper-style figure
show_marker   = true;

%% ===================== Plot style =====================

% Color version: clearer for screen / slide / report.
method_colors = [
    0.00 0.00 0.00;   % POE  black
    0.00 0.25 0.80;   % GPOE blue
    0.80 0.10 0.10;   % MOE  red
    0.10 0.55 0.10;   % BCM  green
    0.55 0.10 0.70    % RBCM purple
];

% If you want black-white paper style, uncomment this:
% method_colors = repmat([0 0 0], numel(Method_list), 1);

method_markers = {'o','^','s','d','v'};
method_lines   = {'-','-','-','-','-'};

line_width  = 1.20;
marker_size = 4.0;

% Do not put marker on every point.
marker_idx = 1:4:numel(M_list);

% Sparse errorbar index, only used when show_errorbar = true.
err_idx = 1:5:numel(M_list);

%% ===================== Read data =====================

data = struct();
MissingCounter = 0;

for ci = 1:numel(Method_list)
    method = Method_list{ci};

    for di = 1:numel(Dataset_list)
        dataset = Dataset_list{di};

        ResultFolder = fullfile(ResultRoot, dataset);

        SMSE_all = nan(numel(M_list), numel(seeds));
        MSLL_all = nan(numel(M_list), numel(seeds));

        for mi = 1:numel(M_list)
            M = M_list(mi);

            for si = 1:numel(seeds)
                seed = seeds(si);

                fname = sprintf('%s_M%d_tr%d_mc%d.mat', method, M, tr_tag, seed);
                fpath = fullfile(ResultFolder, fname);

                if ~exist(fpath, 'file')
                    MissingCounter = MissingCounter + 1;
                    fprintf('[Missing] %s\n', fpath);
                    continue;
                end

                S = load(fpath, 'smse', 'msll');

                if isfield(S, 'smse')
                    SMSE_all(mi, si) = S.smse;
                end

                if isfield(S, 'msll')
                    MSLL_all(mi, si) = S.msll;
                end
            end
        end

        data.(method).(dataset).SMSE_mean = mean(SMSE_all, 2, 'omitnan');
        data.(method).(dataset).SMSE_std  = std(SMSE_all, 0, 2, 'omitnan');

        data.(method).(dataset).MSLL_mean = mean(MSLL_all, 2, 'omitnan');
        data.(method).(dataset).MSLL_std  = std(MSLL_all, 0, 2, 'omitnan');

        data.(method).(dataset).SMSE_all = SMSE_all;
        data.(method).(dataset).MSLL_all = MSLL_all;
    end
end

fprintf('\nData loading finished. Missing files: %d\n', MissingCounter);

%% ===================== Export summary table =====================

rows = {};

for di = 1:numel(Dataset_list)
    dataset = Dataset_list{di};

    for ci = 1:numel(Method_list)
        method = Method_list{ci};
        d = data.(method).(dataset);

        for mi = 1:numel(M_list)
            rows(end+1,:) = { ...
                dataset, ...
                upper(method), ...
                M_list(mi), ...
                d.SMSE_mean(mi), ...
                d.SMSE_std(mi), ...
                d.MSLL_mean(mi), ...
                d.MSLL_std(mi) ...
            }; %#ok<SAGROW>
        end
    end
end

SummaryTable = cell2table(rows, ...
    'VariableNames', {'Dataset','Method','M', ...
    'SMSE_mean','SMSE_std','MSLL_mean','MSLL_std'});

csv_path = fullfile(FigFolder, 'M_ablation_summary_DAC_tr40.csv');
writetable(SummaryTable, csv_path);

fprintf('Summary table saved to:\n%s\n', csv_path);

%% ===================== Print best M summary =====================

fprintf('\n============================================================\n');
fprintf('Best M summary based on minimum SMSE and minimum MSLL\n');
fprintf('============================================================\n');
fprintf('%-12s %-6s | %-22s | %-22s\n', ...
    'Dataset', 'Method', 'Best SMSE', 'Best MSLL');
fprintf('%s\n', repmat('-',1,72));

for di = 1:numel(Dataset_list)
    dataset = Dataset_list{di};

    for ci = 1:numel(Method_list)
        method = Method_list{ci};
        d = data.(method).(dataset);

        [best_smse, idx_smse] = min(d.SMSE_mean);
        [best_msll, idx_msll] = min(d.MSLL_mean);

        fprintf('%-12s %-6s | M=%4d, %.4g       | M=%4d, %.4g\n', ...
            dataset, upper(method), ...
            M_list(idx_smse), best_smse, ...
            M_list(idx_msll), best_msll);
    end
end

fprintf('%s\n', repmat('-',1,72));

%% ===================== Draw figures =====================

draw_tradeoff_figure( ...
    data, Method_list, Method_label, Dataset_list, M_list, ...
    method_colors, method_markers, method_lines, ...
    marker_idx, err_idx, show_errorbar, show_marker, ...
    line_width, marker_size, ...
    [100 2500], ...
    'Effect of the Number of Inducing Points under IP-DAC', ...
    fullfile(FigFolder, 'M_ablation_overall_M100_2500'));

draw_tradeoff_figure( ...
    data, Method_list, Method_label, Dataset_list, M_list, ...
    method_colors, method_markers, method_lines, ...
    marker_idx, err_idx, show_errorbar, show_marker, ...
    line_width, marker_size, ...
    [500 2500], ...
    'Zoomed View: Effect of the Number of Inducing Points under IP-DAC', ...
    fullfile(FigFolder, 'M_ablation_zoom_M500_2500'));

fprintf('\nAll figures saved to:\n%s\n', FigFolder);

%% ========================================================================
%% Local function: draw one 2 x 4 paper-style figure
%% ========================================================================

function draw_tradeoff_figure( ...
    data, Method_list, Method_label, Dataset_list, M_list, ...
    method_colors, method_markers, method_lines, ...
    marker_idx, err_idx, show_errorbar, show_marker, ...
    line_width, marker_size, ...
    x_range, fig_title, save_prefix)

    nD = numel(Dataset_list);
    nM = numel(Method_list);

    x_mask = M_list >= x_range(1) & M_list <= x_range(2);
    x_plot = M_list(x_mask);

    fig = figure('Color','w', ...
        'Units','centimeters', ...
        'Position',[2, 2, 42, 17]);

    tl = tiledlayout(2, nD, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    legend_handles = gobjects(nM, 1);

    for di = 1:nD
        dataset = Dataset_list{di};

        %% -------------------- Top row: SMSE --------------------
        ax1 = nexttile(di);
        hold(ax1, 'on');

        all_smse_vals = [];

        for ci = 1:nM
            method = Method_list{ci};
            d = data.(method).(dataset);

            y = d.SMSE_mean(:);
            e = d.SMSE_std(:);

            x = M_list(:);
            valid = x_mask(:) & isfinite(y) & y > 0;

            xv = x(valid);
            yv = y(valid);
            ev = e(valid);

            if isempty(xv)
                continue;
            end

            h = plot(ax1, xv, yv, ...
                'LineStyle', method_lines{ci}, ...
                'Color', method_colors(ci,:), ...
                'LineWidth', line_width, ...
                'DisplayName', Method_label{ci});

            if show_marker
                local_marker_idx = marker_idx(marker_idx <= numel(xv));
                plot(ax1, xv(local_marker_idx), yv(local_marker_idx), ...
                    'LineStyle','none', ...
                    'Marker', method_markers{ci}, ...
                    'MarkerSize', marker_size, ...
                    'Color', method_colors(ci,:), ...
                    'HandleVisibility','off');
            end

            if show_errorbar
                local_err_idx = err_idx(err_idx <= numel(xv));
                errorbar(ax1, xv(local_err_idx), yv(local_err_idx), ev(local_err_idx), ...
                    'LineStyle','none', ...
                    'Color', method_colors(ci,:), ...
                    'LineWidth', 0.6, ...
                    'CapSize', 2, ...
                    'HandleVisibility','off');
            end

            all_smse_vals = [all_smse_vals; yv]; %#ok<AGROW>

            if di == 1
                legend_handles(ci) = h;
            end
        end

        set(ax1, 'XScale','log', 'YScale','log');
        xlim(ax1, x_range);
        setup_x_ticks_scientific(ax1, x_range);
        xticklabels(ax1, []);

        grid(ax1, 'on');
        box(ax1, 'on');

        title(ax1, dataset, ...
            'FontSize', 11, ...
            'FontWeight','normal', ...
            'Interpreter','none');

        ylabel(ax1, 'SMSE', 'FontSize', 11);

        set(ax1, ...
            'FontSize', 9, ...
            'LineWidth', 0.8, ...
            'TickDir','out', ...
            'TickLabelInterpreter','latex');

        if ~isempty(all_smse_vals)
            all_smse_vals = all_smse_vals(isfinite(all_smse_vals) & all_smse_vals > 0);
            if ~isempty(all_smse_vals)
                ylim(ax1, [min(all_smse_vals)*0.85, max(all_smse_vals)*1.20]);
            end
        end

        %% -------------------- Bottom row: MSLL --------------------
        ax2 = nexttile(nD + di);
        hold(ax2, 'on');

        all_msll_vals = [];

        for ci = 1:nM
            method = Method_list{ci};
            d = data.(method).(dataset);

            y = d.MSLL_mean(:);
            e = d.MSLL_std(:);

            x = M_list(:);
            valid = x_mask(:) & isfinite(y);

            xv = x(valid);
            yv = y(valid);
            ev = e(valid);

            if isempty(xv)
                continue;
            end

            plot(ax2, xv, yv, ...
                'LineStyle', method_lines{ci}, ...
                'Color', method_colors(ci,:), ...
                'LineWidth', line_width, ...
                'HandleVisibility','off');

            if show_marker
                local_marker_idx = marker_idx(marker_idx <= numel(xv));
                plot(ax2, xv(local_marker_idx), yv(local_marker_idx), ...
                    'LineStyle','none', ...
                    'Marker', method_markers{ci}, ...
                    'MarkerSize', marker_size, ...
                    'Color', method_colors(ci,:), ...
                    'HandleVisibility','off');
            end

            if show_errorbar
                local_err_idx = err_idx(err_idx <= numel(xv));
                errorbar(ax2, xv(local_err_idx), yv(local_err_idx), ev(local_err_idx), ...
                    'LineStyle','none', ...
                    'Color', method_colors(ci,:), ...
                    'LineWidth', 0.6, ...
                    'CapSize', 2, ...
                    'HandleVisibility','off');
            end

            all_msll_vals = [all_msll_vals; yv]; %#ok<AGROW>
        end

        set(ax2, 'XScale','log');
        xlim(ax2, x_range);
        setup_x_ticks_scientific(ax2, x_range);

        grid(ax2, 'on');
        box(ax2, 'on');

        xlabel(ax2, '$M$', ...
            'FontSize', 11, ...
            'Interpreter','latex');

        ylabel(ax2, 'MSLL', 'FontSize', 11);

        set(ax2, ...
            'FontSize', 9, ...
            'LineWidth', 0.8, ...
            'TickDir','out', ...
            'TickLabelInterpreter','latex');

        if ~isempty(all_msll_vals)
            all_msll_vals = all_msll_vals(isfinite(all_msll_vals));
            if ~isempty(all_msll_vals)
                pad = 0.08 * (max(all_msll_vals) - min(all_msll_vals) + eps);
                ylim(ax2, [min(all_msll_vals)-pad, max(all_msll_vals)+pad]);
            end
        end
    end

    lgd = legend(legend_handles, Method_label, ...
        'Orientation','horizontal', ...
        'NumColumns', numel(Method_list), ...
        'FontSize', 10, ...
        'Box','off');

    lgd.Layout.Tile = 'south';

    title(tl, fig_title, ...
        'FontSize', 13, ...
        'FontWeight','normal');

    exportgraphics(fig, [save_prefix '.png'], 'Resolution', 300);
    exportgraphics(fig, [save_prefix '.pdf'], 'ContentType','vector');
    savefig(fig, [save_prefix '.fig']);
end

%% ========================================================================
%% Local function: compact scientific x tick labels
%% ========================================================================

function setup_x_ticks_scientific(ax, x_range)

    if x_range(1) <= 100
        xticks(ax, [100 500 1000 2500]);
        xticklabels(ax, {'$10^2$', '$5{\times}10^2$', '$10^3$', '$2.5{\times}10^3$'});
    else
        xticks(ax, [500 1000 1500 2000 2500]);
        xticklabels(ax, {'$5{\times}10^2$', '$10^3$', '$1.5{\times}10^3$', ...
                         '$2{\times}10^3$', '$2.5{\times}10^3$'});
    end

    ax.XMinorTick = 'off';
end