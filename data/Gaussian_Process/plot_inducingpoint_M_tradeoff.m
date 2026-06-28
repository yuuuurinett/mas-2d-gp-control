%% plot_inducingpoint_M_tradeoff.m
% Paper-style M-ablation summary:
%   Fixed setting: IP-DAC / event-triggered consensus
%   M = 100:100:2500
%   Metrics: SMSE, MSLL, and train time
%   Mean/std over Monte Carlo seeds
%
% Outputs:
%   1) Full M-ablation summary CSV
%   2) Low / Medium / High M tables
%   3) Selected inducing-point number table
%   4) Paper-style figures
%
% Notes:
%   Low    = M=500
%   Medium = M=1500
%   High   = M=2500

clear; clc; close all;

%% ===================== User settings =====================

Method_list  = {'poe','gpoe','moe','bcm','rbcm'};
Method_label = {'POE','GPOE','MOE','BCM','RBCM'};

Dataset_list = {'KIN40K','POL','PUMADYN32NM','SARCOS'};

M_list      = 100:100:2500;
seeds       = 1:3;
train_ratio = 0.4;
tr_tag      = round(train_ratio * 100);

ProjectRoot = fileparts(mfilename('fullpath'));

% If needed, force it manually:
% ProjectRoot = 'C:\Users\Yurou Du\Desktop\mas-2d-gp-control\data\Gaussian_Process';

ResultRoot = fullfile(ProjectRoot, 'Result', 'Dataset');

FigFolder = fullfile(ProjectRoot, 'Result', 'Figures', 'M_ablation_pretty');
if ~exist(FigFolder, 'dir')
    mkdir(FigFolder);
end

show_errorbar = false;
show_marker   = true;

% M is displayed in units of 10^3 on the x-axis.
x_unit = 1000;

% Overall and zoom ranges in original M units.
overall_range = [100 2500];
zoom_range    = [1000 2500];

%% ===================== Plot style =====================

method_colors = [
    0.00 0.00 0.00;   % POE  black
    0.00 0.25 0.80;   % GPOE blue
    0.80 0.10 0.10;   % MOE  red
    0.10 0.55 0.10;   % BCM  green
    0.55 0.10 0.70    % RBCM purple
];

method_markers = {'o','^','s','d','v'};
method_lines   = {'-','-','-','-','-'};

line_width  = 1.35;
marker_size = 4.2;

marker_idx = 1:4:numel(M_list);
err_idx    = 1:5:numel(M_list);

%% ===================== Read data =====================

data = struct();
MissingCounter = 0;

for ci = 1:numel(Method_list)
    method = Method_list{ci};

    for di = 1:numel(Dataset_list)
        dataset = Dataset_list{di};

        ResultFolder = fullfile(ResultRoot, dataset);

        SMSE_all      = nan(numel(M_list), numel(seeds));
        MSLL_all      = nan(numel(M_list), numel(seeds));
        TrainTime_all = nan(numel(M_list), numel(seeds));   % ms/pt

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

                S = load(fpath, 'smse', 'msll', ...
                    't_train_per_point', 't_train_total', 'N_train');

                if isfield(S, 'smse')
                    SMSE_all(mi, si) = S.smse;
                end

                if isfield(S, 'msll')
                    MSLL_all(mi, si) = S.msll;
                end

                % Training time in ms/pt.
                % New result files should contain t_train_per_point directly.
                % For older result files, fall back to t_train_total / N_train.
                if isfield(S, 't_train_per_point')
                    TrainTime_all(mi, si) = S.t_train_per_point;
                elseif isfield(S, 't_train_total') && isfield(S, 'N_train') && S.N_train > 0
                    TrainTime_all(mi, si) = (S.t_train_total / S.N_train) * 1000;
                end
            end
        end

        data.(method).(dataset).SMSE_mean = mean(SMSE_all, 2, 'omitnan');
        data.(method).(dataset).SMSE_std  = std(SMSE_all, 0, 2, 'omitnan');

        data.(method).(dataset).MSLL_mean = mean(MSLL_all, 2, 'omitnan');
        data.(method).(dataset).MSLL_std  = std(MSLL_all, 0, 2, 'omitnan');

        data.(method).(dataset).TrainTime_mean = mean(TrainTime_all, 2, 'omitnan');
        data.(method).(dataset).TrainTime_std  = std(TrainTime_all, 0, 2, 'omitnan');

        data.(method).(dataset).SMSE_all      = SMSE_all;
        data.(method).(dataset).MSLL_all      = MSLL_all;
        data.(method).(dataset).TrainTime_all = TrainTime_all;
    end
end

fprintf('\nData loading finished. Missing files: %d\n', MissingCounter);

%% ===================== Export full summary table =====================

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
                d.MSLL_std(mi), ...
                d.TrainTime_mean(mi), ...
                d.TrainTime_std(mi) ...
            }; %#ok<SAGROW>
        end
    end
end

SummaryTable = cell2table(rows, ...
    'VariableNames', {'Dataset','Method','M', ...
    'SMSE_mean','SMSE_std','MSLL_mean','MSLL_std', ...
    'TrainTime_mspt_mean','TrainTime_mspt_std'});

csv_path = fullfile(FigFolder, 'M_ablation_summary_DAC_tr40.csv');
writetable(SummaryTable, csv_path);
fprintf('Summary table saved to:\n%s\n', csv_path);

%% ===================== Print best M summary =====================

fprintf('\n============================================================\n');
fprintf('Best M summary based on minimum SMSE, minimum MSLL, and train time\n');
fprintf('============================================================\n');

fprintf('%-12s %-6s | %-22s | %-22s | %-18s\n', ...
    'Dataset', 'Method', 'Best SMSE', 'Best MSLL', 'Time@BestMSLL');
fprintf('%s\n', repmat('-',1,96));

for di = 1:numel(Dataset_list)
    dataset = Dataset_list{di};

    for ci = 1:numel(Method_list)
        method = Method_list{ci};
        d = data.(method).(dataset);

        [best_smse, idx_smse] = min(d.SMSE_mean);
        [best_msll, idx_msll] = min(d.MSLL_mean);

        time_best_msll = d.TrainTime_mean(idx_msll);

        fprintf('%-12s %-6s | M=%4d, %.4g       | M=%4d, %.4g       | %.4g ms/pt\n', ...
            dataset, upper(method), ...
            M_list(idx_smse), best_smse, ...
            M_list(idx_msll), best_msll, ...
            time_best_msll);
    end
end

fprintf('%s\n', repmat('-',1,72));

%% ===================== Low / Medium / High and selected M tables =====================

M_level_name = {'Low','Medium','High'};
M_level      = [500 1500 2500];

fprintf('\n\n============================================================\n');
fprintf('Low / Medium / High M trade-off tables\n');
fprintf('Averaged over seeds = 1:3\n');
fprintf('============================================================\n\n');

%% ------------------------------------------------------------------------
% Detailed table sorted as Dataset -> Level/M -> Method
% -------------------------------------------------------------------------

LevelRows = {};

fprintf('## Detailed Low / Medium / High Table\n\n');
fprintf('Averaged over seeds = 1:3. Sorted as Dataset -> Level/M -> Method.\n\n');
fprintf('| Dataset | Level | M | Method | Train Time (ms/pt) | SMSE | MSLL |\n');
fprintf('|---|---|---:|---|---:|---:|---:|\n');

for di = 1:numel(Dataset_list)
    dataset = Dataset_list{di};

    for li = 1:numel(M_level)
        level_name = M_level_name{li};
        M = M_level(li);
        idx_M = find(M_list == M, 1);

        if isempty(idx_M)
            warning('M=%d is not found in M_list.', M);
            continue;
        end

        for ci = 1:numel(Method_list)
            method = Method_list{ci};
            method_label = Method_label{ci};
            d = data.(method).(dataset);

            train_time_mean = d.TrainTime_mean(idx_M);
            train_time_std  = d.TrainTime_std(idx_M);

            smse_mean = d.SMSE_mean(idx_M);
            smse_std  = d.SMSE_std(idx_M);

            msll_mean = d.MSLL_mean(idx_M);
            msll_std  = d.MSLL_std(idx_M);

            fprintf('| %s | %s | %d | %s | %.4g | %.4g | %.4g |\n', ...
                dataset, level_name, M, method_label, ...
                train_time_mean, smse_mean, msll_mean);

            LevelRows(end+1,:) = { ...
                dataset, level_name, M, method_label, ...
                train_time_mean, train_time_std, ...
                smse_mean, smse_std, ...
                msll_mean, msll_std ...
            }; %#ok<SAGROW>
        end
    end
end

LevelTable = cell2table(LevelRows, ...
    'VariableNames', {'Dataset','Level','M','Method', ...
    'TrainTime_mean_mspt','TrainTime_std_mspt', ...
    'SMSE_mean','SMSE_std', ...
    'MSLL_mean','MSLL_std'});

level_csv_path = fullfile(FigFolder, 'M_level_tradeoff_table_dataset_level_method.csv');
writetable(LevelTable, level_csv_path);

fprintf('\nDetailed CSV table saved to:\n%s\n', level_csv_path);

%% ------------------------------------------------------------------------
% Compact table averaged over methods
% -------------------------------------------------------------------------

fprintf('\n\n## Compact Low / Medium / High Table Averaged over Methods\n\n');
fprintf('This compact table averages POE, GPOE, MOE, BCM, and RBCM for each dataset and M level.\n\n');
fprintf('| Dataset | Level | M | Avg Train Time (ms/pt) | Avg SMSE | Avg MSLL | Observation |\n');
fprintf('|---|---|---:|---:|---:|---:|---|\n');

CompactRows = {};

for di = 1:numel(Dataset_list)
    dataset = Dataset_list{di};

    for li = 1:numel(M_level)
        level_name = M_level_name{li};
        M = M_level(li);
        idx_M = find(M_list == M, 1);

        if isempty(idx_M)
            warning('M=%d is not found in M_list.', M);
            continue;
        end

        time_vals = [];
        smse_vals = [];
        msll_vals = [];

        for ci = 1:numel(Method_list)
            method = Method_list{ci};
            d = data.(method).(dataset);

            time_vals(end+1) = d.TrainTime_mean(idx_M); %#ok<SAGROW>
            smse_vals(end+1) = d.SMSE_mean(idx_M); %#ok<SAGROW>
            msll_vals(end+1) = d.MSLL_mean(idx_M); %#ok<SAGROW>
        end

        avg_time = mean(time_vals, 'omitnan');
        avg_smse = mean(smse_vals, 'omitnan');
        avg_msll = mean(msll_vals, 'omitnan');

        observation = make_level_observation(dataset, level_name);

        fprintf('| %s | %s | %d | %.4g | %.4g | %.4g | %s |\n', ...
            dataset, level_name, M, avg_time, avg_smse, avg_msll, observation);

        CompactRows(end+1,:) = { ...
            dataset, level_name, M, ...
            avg_time, avg_smse, avg_msll, observation ...
        }; %#ok<SAGROW>
    end
end

CompactTable = cell2table(CompactRows, ...
    'VariableNames', {'Dataset','Level','M', ...
    'AvgTrainTime_mspt','AvgSMSE','AvgMSLL','Observation'});

compact_csv_path = fullfile(FigFolder, 'M_level_tradeoff_compact_avg_methods.csv');
writetable(CompactTable, compact_csv_path);

fprintf('\nCompact CSV table saved to:\n%s\n', compact_csv_path);

%% ------------------------------------------------------------------------
% Selected M decision table
% -------------------------------------------------------------------------

SelectedRows = {
    'KIN40K',      1800, 'Balanced choice: best MSLL region with strong SMSE improvement; M=2500 improves SMSE but worsens MSLL and costs more.';
    'POL',         2000, 'Balanced choice: near-best SMSE and best MSLL region; M=2500 gives slightly lower SMSE but weaker MSLL.';
    'PUMADYN32NM', 1500, 'Saturated dataset: M=1500 is sufficient; larger M increases training time with negligible metric gain.';
    'SARCOS',      2500, 'Largest tested M is justified: both SMSE and MSLL continue improving up to M=2500.';
};

SelectedTable = cell2table(SelectedRows, ...
    'VariableNames', {'Dataset','SelectedM','Reason'});

fprintf('\n\n## Selected Number of Inducing Points\n\n');
fprintf('| Dataset | Selected M | Reason |\n');
fprintf('|---|---:|---|\n');

for r = 1:height(SelectedTable)
    fprintf('| %s | %d | %s |\n', ...
        SelectedTable.Dataset{r}, ...
        SelectedTable.SelectedM(r), ...
        SelectedTable.Reason{r});
end

selected_csv_path = fullfile(FigFolder, 'M_selected_inducing_points.csv');
writetable(SelectedTable, selected_csv_path);

fprintf('\nSelected M CSV table saved to:\n%s\n', selected_csv_path);

%% ------------------------------------------------------------------------
% Save all tables to Markdown file
% -------------------------------------------------------------------------

md_path = fullfile(FigFolder, 'M_level_tradeoff_and_selection_tables.md');
fid = fopen(md_path, 'w');

fprintf(fid, '# M-level Trade-off and Inducing Point Selection\n\n');
fprintf(fid, 'Averaged over seeds = 1:3. Low = M=500, Medium = M=1500, High = M=2500.\n\n');

fprintf(fid, '## Recommended Number of Inducing Points\n\n');
fprintf(fid, '| Dataset | Selected M | Reason |\n');
fprintf(fid, '|---|---:|---|\n');

for r = 1:height(SelectedTable)
    fprintf(fid, '| %s | %d | %s |\n', ...
        SelectedTable.Dataset{r}, ...
        SelectedTable.SelectedM(r), ...
        SelectedTable.Reason{r});
end

fprintf(fid, '\n\n## Compact Low / Medium / High Table Averaged over Methods\n\n');
fprintf(fid, '| Dataset | Level | M | Avg Train Time (ms/pt) | Avg SMSE | Avg MSLL | Observation |\n');
fprintf(fid, '|---|---|---:|---:|---:|---:|---|\n');

for r = 1:height(CompactTable)
    fprintf(fid, '| %s | %s | %d | %.4g | %.4g | %.4g | %s |\n', ...
        CompactTable.Dataset{r}, ...
        CompactTable.Level{r}, ...
        CompactTable.M(r), ...
        CompactTable.AvgTrainTime_mspt(r), ...
        CompactTable.AvgSMSE(r), ...
        CompactTable.AvgMSLL(r), ...
        CompactTable.Observation{r});
end

fprintf(fid, '\n\n## Detailed Grouped Table by Dataset, M Level, and Method\n\n');
fprintf(fid, 'This detailed table is sorted as Dataset -> Level/M -> Method.\n\n');

fprintf(fid, '<table>\n');

fprintf(fid, '<tr>\n');
fprintf(fid, '<th>Dataset</th>\n');
fprintf(fid, '<th>Level</th>\n');
fprintf(fid, '<th>M</th>\n');
fprintf(fid, '<th>Method</th>\n');
fprintf(fid, '<th>Train Time (ms/pt)</th>\n');
fprintf(fid, '<th>SMSE</th>\n');
fprintf(fid, '<th>MSLL</th>\n');
fprintf(fid, '</tr>\n');

for di = 1:numel(Dataset_list)
    dataset = Dataset_list{di};

    for li = 1:numel(M_level)
        level_name = M_level_name{li};
        M = M_level(li);
        idx_M = find(M_list == M, 1);

        if isempty(idx_M)
            continue;
        end

        for ci = 1:numel(Method_list)
            method = Method_list{ci};
            method_label = Method_label{ci};
            d = data.(method).(dataset);

            train_time_mean = d.TrainTime_mean(idx_M);
            smse_mean = d.SMSE_mean(idx_M);
            msll_mean = d.MSLL_mean(idx_M);

            fprintf(fid, '<tr>\n');

            % Dataset spans 3 levels x 5 methods = 15 rows.
            if li == 1 && ci == 1
                fprintf(fid, '<td rowspan="15"><b>%s</b></td>\n', dataset);
            end

            % Level and M each span 5 method rows.
            if ci == 1
                fprintf(fid, '<td rowspan="5"><b>%s</b></td>\n', level_name);
                fprintf(fid, '<td rowspan="5">%d</td>\n', M);
            end

            fprintf(fid, '<td>%s</td>\n', method_label);
            fprintf(fid, '<td>%.4g</td>\n', train_time_mean);
            fprintf(fid, '<td>%.4g</td>\n', smse_mean);
            fprintf(fid, '<td>%.4g</td>\n', msll_mean);

            fprintf(fid, '</tr>\n');
        end
    end
end

fprintf(fid, '</table>\n');

fclose(fid);

fprintf('\nMarkdown tables saved to:\n%s\n', md_path);

%% ===================== Draw figures =====================

draw_tradeoff_figure( ...
    data, Method_list, Method_label, Dataset_list, M_list, x_unit, ...
    method_colors, method_markers, method_lines, ...
    marker_idx, err_idx, show_errorbar, show_marker, ...
    line_width, marker_size, ...
    overall_range, false, ...
    'Effect of the Number of Inducing Points under IP-DAC', ...
    fullfile(FigFolder, 'M_ablation_overall_M100_2500_pretty'));

draw_tradeoff_figure( ...
    data, Method_list, Method_label, Dataset_list, M_list, x_unit, ...
    method_colors, method_markers, method_lines, ...
    marker_idx, err_idx, show_errorbar, show_marker, ...
    line_width, marker_size, ...
    zoom_range, true, ...
    'Zoomed View for $M \geq 10^3$ under IP-DAC', ...
    fullfile(FigFolder, 'M_ablation_zoom_M1000_2500_pretty'));

draw_time_tradeoff_figure( ...
    data, Method_list, Method_label, Dataset_list, M_list, ...
    method_colors, method_markers, method_lines, ...
    line_width, marker_size, ...
    overall_range, false, ...
    'Train-Time Trade-off under IP-DAC', ...
    fullfile(FigFolder, 'M_ablation_train_time_tradeoff_overall'));

draw_time_tradeoff_figure( ...
    data, Method_list, Method_label, Dataset_list, M_list, ...
    method_colors, method_markers, method_lines, ...
    line_width, marker_size, ...
    zoom_range, true, ...
    'Zoomed Train-Time Trade-off for $M \geq 10^3$ under IP-DAC', ...
    fullfile(FigFolder, 'M_ablation_train_time_tradeoff_zoom_M1000_2500'));

fprintf('\nAll figures saved to:\n%s\n', FigFolder);

%% ========================================================================
%% Local function: draw one 2 x 4 paper-style figure
%% ========================================================================

function draw_tradeoff_figure( ...
    data, Method_list, Method_label, Dataset_list, M_list, x_unit, ...
    method_colors, method_markers, method_lines, ...
    marker_idx, err_idx, show_errorbar, show_marker, ...
    line_width, marker_size, ...
    x_range_M, is_zoom, fig_title, save_prefix)

    nD = numel(Dataset_list);
    nM = numel(Method_list);

    x_mask = M_list >= x_range_M(1) & M_list <= x_range_M(2);
    x_plot_all = M_list(:) / x_unit;
    x_range = x_range_M / x_unit;

    if is_zoom
        fig_pos = [2, 2, 42, 17.5];
        smse_min_range = 0.006;
        msll_min_range = 0.040;
    else
        fig_pos = [2, 2, 42, 17.5];
        smse_min_range = 0.015;
        msll_min_range = 0.080;
    end

    fig = figure('Color','w', ...
        'Units','centimeters', ...
        'Position',fig_pos);

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

            valid = x_mask(:) & isfinite(y) & y > 0;

            xv = x_plot_all(valid);
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
                local_marker_idx = choose_marker_indices(numel(xv), is_zoom);
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

        set(ax1, 'YScale','log');
        xlim(ax1, x_range);
        setup_x_ticks_scaled(ax1, x_range_M, x_unit, false);
        xticklabels(ax1, []);

        grid(ax1, 'on');
        box(ax1, 'on');

        title(ax1, dataset, ...
            'FontSize', 12, ...
            'FontWeight','normal', ...
            'Interpreter','none');

        ylabel(ax1, 'SMSE', 'FontSize', 12);

        set(ax1, ...
            'FontSize', 10, ...
            'LineWidth', 0.85, ...
            'TickDir','out', ...
            'TickLabelInterpreter','latex');

        apply_log_ylim(ax1, all_smse_vals, 0.07, smse_min_range);

        %% -------------------- Bottom row: MSLL --------------------
        ax2 = nexttile(nD + di);
        hold(ax2, 'on');

        all_msll_vals = [];

        for ci = 1:nM
            method = Method_list{ci};
            d = data.(method).(dataset);

            y = d.MSLL_mean(:);
            e = d.MSLL_std(:);

            valid = x_mask(:) & isfinite(y);

            xv = x_plot_all(valid);
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
                local_marker_idx = choose_marker_indices(numel(xv), is_zoom);
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

        xlim(ax2, x_range);
        setup_x_ticks_scaled(ax2, x_range_M, x_unit, true);

        grid(ax2, 'on');
        box(ax2, 'on');

        xlabel(ax2, '$M\;(\times 10^3)$', ...
            'FontSize', 12, ...
            'Interpreter','latex');

        ylabel(ax2, 'MSLL', 'FontSize', 12);

        set(ax2, ...
            'FontSize', 10, ...
            'LineWidth', 0.85, ...
            'TickDir','out', ...
            'TickLabelInterpreter','latex');

        apply_linear_ylim(ax2, all_msll_vals, 0.08, msll_min_range);
    end

    lgd = legend(legend_handles, Method_label, ...
        'Orientation','horizontal', ...
        'NumColumns', numel(Method_list), ...
        'FontSize', 11, ...
        'Box','off');

    lgd.Layout.Tile = 'south';

    title(tl, fig_title, ...
        'FontSize', 15, ...
        'FontWeight','normal', ...
        'Interpreter','latex');

    exportgraphics(fig, [save_prefix '.png'], 'Resolution', 300);
    exportgraphics(fig, [save_prefix '.pdf'], 'ContentType','vector');
    savefig(fig, [save_prefix '.fig']);
end

%% ========================================================================
%% Local function: draw train-time trade-off figure
%% ========================================================================

function draw_time_tradeoff_figure( ...
    data, Method_list, Method_label, Dataset_list, M_list, ...
    method_colors, method_markers, method_lines, ...
    line_width, marker_size, ...
    x_range_M, is_zoom, fig_title, save_prefix)

    nD = numel(Dataset_list);
    nM = numel(Method_list);

    x_mask = M_list >= x_range_M(1) & M_list <= x_range_M(2);

    if is_zoom
        fig_pos = [2, 2, 42, 17.5];
        smse_min_range = 0.006;
        msll_min_range = 0.040;
    else
        fig_pos = [2, 2, 42, 17.5];
        smse_min_range = 0.015;
        msll_min_range = 0.080;
    end

    fig = figure('Color','w', ...
        'Units','centimeters', ...
        'Position',fig_pos);

    tl = tiledlayout(2, nD, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    legend_handles = gobjects(nM, 1);

    for di = 1:nD
        dataset = Dataset_list{di};

        %% -------------------- Top row: Train time vs SMSE --------------------
        ax1 = nexttile(di);
        hold(ax1, 'on');

        all_smse_vals = [];
        all_time_vals = [];

        for ci = 1:nM
            method = Method_list{ci};
            d = data.(method).(dataset);

            x = d.TrainTime_mean(:);   % ms/pt
            y = d.SMSE_mean(:);

            valid = x_mask(:) & isfinite(x) & isfinite(y) & x > 0 & y > 0;
            xv = x(valid);
            yv = y(valid);
            Mv = M_list(valid);

            if isempty(xv)
                continue;
            end

            [~, order] = sort(Mv);
            xv = xv(order);
            yv = yv(order);

            h = plot(ax1, xv, yv, ...
                'LineStyle', method_lines{ci}, ...
                'Color', method_colors(ci,:), ...
                'LineWidth', line_width, ...
                'DisplayName', Method_label{ci});

            local_marker_idx = choose_marker_indices(numel(xv), is_zoom);
            plot(ax1, xv(local_marker_idx), yv(local_marker_idx), ...
                'LineStyle','none', ...
                'Marker', method_markers{ci}, ...
                'MarkerSize', marker_size, ...
                'Color', method_colors(ci,:), ...
                'HandleVisibility','off');

            all_smse_vals = [all_smse_vals; yv]; %#ok<AGROW>
            all_time_vals = [all_time_vals; xv]; %#ok<AGROW>

            if di == 1
                legend_handles(ci) = h;
            end
        end

        set(ax1, 'YScale','log');
        apply_linear_xlim(ax1, all_time_vals, 0.08, 0.02);
        apply_log_ylim(ax1, all_smse_vals, 0.07, smse_min_range);

        grid(ax1, 'on');
        box(ax1, 'on');

        title(ax1, dataset, ...
            'FontSize', 12, ...
            'FontWeight','normal', ...
            'Interpreter','none');

        ylabel(ax1, 'SMSE', 'FontSize', 12);
        xticklabels(ax1, []);

        set(ax1, ...
            'FontSize', 10, ...
            'LineWidth', 0.85, ...
            'TickDir','out', ...
            'TickLabelInterpreter','latex');

        %% -------------------- Bottom row: Train time vs MSLL --------------------
        ax2 = nexttile(nD + di);
        hold(ax2, 'on');

        all_msll_vals = [];
        all_time_vals = [];

        for ci = 1:nM
            method = Method_list{ci};
            d = data.(method).(dataset);

            x = d.TrainTime_mean(:);   % ms/pt
            y = d.MSLL_mean(:);

            valid = x_mask(:) & isfinite(x) & isfinite(y) & x > 0;
            xv = x(valid);
            yv = y(valid);
            Mv = M_list(valid);

            if isempty(xv)
                continue;
            end

            [~, order] = sort(Mv);
            xv = xv(order);
            yv = yv(order);

            plot(ax2, xv, yv, ...
                'LineStyle', method_lines{ci}, ...
                'Color', method_colors(ci,:), ...
                'LineWidth', line_width, ...
                'HandleVisibility','off');

            local_marker_idx = choose_marker_indices(numel(xv), is_zoom);
            plot(ax2, xv(local_marker_idx), yv(local_marker_idx), ...
                'LineStyle','none', ...
                'Marker', method_markers{ci}, ...
                'MarkerSize', marker_size, ...
                'Color', method_colors(ci,:), ...
                'HandleVisibility','off');

            all_msll_vals = [all_msll_vals; yv]; %#ok<AGROW>
            all_time_vals = [all_time_vals; xv]; %#ok<AGROW>
        end

        apply_linear_xlim(ax2, all_time_vals, 0.08, 0.02);
        apply_linear_ylim(ax2, all_msll_vals, 0.08, msll_min_range);

        grid(ax2, 'on');
        box(ax2, 'on');

        xlabel(ax2, 'Train Time (ms/pt)', ...
            'FontSize', 12, ...
            'Interpreter','none');

        ylabel(ax2, 'MSLL', 'FontSize', 12);

        set(ax2, ...
            'FontSize', 10, ...
            'LineWidth', 0.85, ...
            'TickDir','out', ...
            'TickLabelInterpreter','latex');
    end

    lgd = legend(legend_handles, Method_label, ...
        'Orientation','horizontal', ...
        'NumColumns', numel(Method_list), ...
        'FontSize', 11, ...
        'Box','off');

    lgd.Layout.Tile = 'south';

    title(tl, fig_title, ...
        'FontSize', 15, ...
        'FontWeight','normal', ...
        'Interpreter','latex');

    exportgraphics(fig, [save_prefix '.png'], 'Resolution', 300);
    exportgraphics(fig, [save_prefix '.pdf'], 'ContentType','vector');
    savefig(fig, [save_prefix '.fig']);
end

%% ========================================================================
%% Local function: simple x tick labels in units of 10^3
%% ========================================================================

function setup_x_ticks_scaled(ax, x_range_M, x_unit, show_labels)

    if x_range_M(1) <= 100
        ticks_M = [100 500 1000 1500 2000 2500];
        labels  = {'0.1','0.5','1','1.5','2','2.5'};
    elseif x_range_M(1) <= 500
        ticks_M = [500 1000 1500 2000 2500];
        labels  = {'0.5','1','1.5','2','2.5'};
    else
        ticks_M = [1000 1500 2000 2500];
        labels  = {'1','1.5','2','2.5'};
    end

    xticks(ax, ticks_M / x_unit);

    if show_labels
        xticklabels(ax, labels);
    else
        xticklabels(ax, []);
    end

    ax.XMinorTick = 'off';
end

%% ========================================================================
%% Local function: axis helpers
%% ========================================================================

function apply_linear_xlim(ax, x_values, pad_ratio, min_range)
    x_values = x_values(:);
    x_values = x_values(isfinite(x_values));

    if isempty(x_values)
        return;
    end

    x_min = min(x_values);
    x_max = max(x_values);
    x_range = x_max - x_min;

    if x_range < min_range
        x_center = 0.5 * (x_min + x_max);
        xlim(ax, [max(0, x_center - min_range/2), x_center + min_range/2]);
    else
        pad = pad_ratio * x_range;
        xlim(ax, [max(0, x_min - pad), x_max + pad]);
    end
end

function apply_log_ylim(ax, y_values, log_pad, min_log_span)
    y_values = y_values(:);
    y_values = y_values(isfinite(y_values) & y_values > 0);

    if isempty(y_values)
        return;
    end

    log_y = log10(y_values);
    lo = min(log_y);
    hi = max(log_y);
    span = hi - lo;

    if span < min_log_span
        center = 0.5 * (lo + hi);
        lo = center - min_log_span/2;
        hi = center + min_log_span/2;
    else
        lo = lo - log_pad * span;
        hi = hi + log_pad * span;
    end

    ylim(ax, 10.^[lo hi]);
end

function apply_linear_ylim(ax, y_values, pad_ratio, min_range)
    y_values = y_values(:);
    y_values = y_values(isfinite(y_values));

    if isempty(y_values)
        return;
    end

    y_min = min(y_values);
    y_max = max(y_values);
    y_range = y_max - y_min;

    if y_range < min_range
        y_center = 0.5 * (y_min + y_max);
        ylim(ax, [y_center - min_range/2, y_center + min_range/2]);
    else
        pad = pad_ratio * y_range;
        ylim(ax, [y_min - pad, y_max + pad]);
    end
end

function idx = choose_marker_indices(n, is_zoom)
    if n <= 0
        idx = [];
        return;
    end

    if is_zoom
        idx = unique([1, round(linspace(1,n,min(5,n))), n]);
    else
        idx = unique([1, round(linspace(1,n,min(6,n))), n]);
    end
end

function observation = make_level_observation(dataset, level_name)
    switch dataset
        case 'KIN40K'
            switch level_name
                case 'Low'
                    observation = 'Underfitting';
                case 'Medium'
                    observation = 'Good trade-off';
                case 'High'
                    observation = 'Better SMSE, weaker MSLL';
                otherwise
                    observation = '';
            end

        case 'POL'
            switch level_name
                case 'Low'
                    observation = 'Underfitting';
                case 'Medium'
                    observation = 'Good trade-off';
                case 'High'
                    observation = 'Better SMSE, slightly weaker MSLL';
                otherwise
                    observation = '';
            end

        case 'PUMADYN32NM'
            switch level_name
                case 'Low'
                    observation = 'Already saturated';
                case 'Medium'
                    observation = 'No meaningful gain';
                case 'High'
                    observation = 'Extra cost only';
                otherwise
                    observation = '';
            end

        case 'SARCOS'
            switch level_name
                case 'Low'
                    observation = 'Good but not saturated';
                case 'Medium'
                    observation = 'Better';
                case 'High'
                    observation = 'Best overall';
                otherwise
                    observation = '';
            end

        otherwise
            observation = '';
    end
end