%% export_summary_tables_to_html.m
% Export editable Word-friendly HTML tables from All_Methods_Summary.mat.
%
% This is NOT image export.
% The generated .html file contains real HTML tables.
% You can:
%   1) open the HTML file in a browser and copy/paste into Word, or
%   2) open the HTML file directly with Microsoft Word.
%
% Layout:
%   For each metric:
%       For each aggregation rule:
%           rows    = LoG / CEN / IP-DAC / IP-AC / TP-DAC / TP-AC / NBR
%           columns = KIN40K / POL / PUMA / SARCOS
%
% Best value in each dataset column is highlighted.
% Lower is better for all current metrics.

clear; clc;

ProjectRoot = fileparts(mfilename('fullpath'));
cd(ProjectRoot);

SummaryPath = fullfile(ProjectRoot, 'Result', 'Dataset', 'All_Methods_Summary.mat');

if ~isfile(SummaryPath)
    error('Summary file not found:\n%s\nRun run_mc_dataset first.', SummaryPath);
end

S = load(SummaryPath, ...
    'mean_results', 'std_results', ...
    'datasets', 'methods_names', 'metrics_names', ...
    'train_ratio', 'n_mc');

mean_results  = S.mean_results;
std_results   = S.std_results;
datasets      = S.datasets;
methods_names = S.methods_names;
metrics_names = S.metrics_names;
train_ratio   = S.train_ratio;
n_mc          = S.n_mc;

agg_list     = {'MOE','GPOE','POE','BCM','RBCM'};
row_families = {'LoG','CEN','IP-DAC','IP-AC','TP-DAC','TP-AC','NBR'};

dataset_labels = datasets;
dataset_labels(strcmp(dataset_labels, 'PUMADYN32NM')) = {'PUMA'};

OutFolder = fullfile(ProjectRoot, 'Result', 'Dataset');
if ~exist(OutFolder, 'dir')
    mkdir(OutFolder);
end

HtmlPath = fullfile(OutFolder, 'All_Methods_Summary_editable_tables.html');

fid = fopen(HtmlPath, 'w');
if fid < 0
    error('Cannot open HTML file for writing:\n%s', HtmlPath);
end

fprintf(fid, '<!DOCTYPE html>\n');
fprintf(fid, '<html>\n<head>\n<meta charset="UTF-8">\n');
fprintf(fid, '<title>All Methods Summary</title>\n');

% CSS for Word-friendly editable tables.
fprintf(fid, '<style>\n');
fprintf(fid, 'body { font-family: "Times New Roman", serif; margin: 24px; }\n');
fprintf(fid, 'h1 { color: #b00000; font-size: 22pt; text-align: center; margin-top: 28px; }\n');
fprintf(fid, 'h2 { color: #b00000; font-size: 18pt; text-align: center; margin-top: 24px; margin-bottom: 8px; }\n');
fprintf(fid, 'p { font-size: 11pt; }\n');
fprintf(fid, 'table { border-collapse: collapse; margin: 10px auto 24px auto; width: 92%%; table-layout: fixed; }\n');
fprintf(fid, 'th, td { border: 1px solid black; padding: 5px 7px; text-align: center; font-size: 10.5pt; }\n');
fprintf(fid, 'th { background-color: #cfe8f3; font-weight: bold; }\n');
fprintf(fid, '.agg { background-color: #cfe8f3; font-weight: bold; font-size: 13pt; }\n');
fprintf(fid, '.rowhead { font-weight: bold; background-color: #f7f7f7; }\n');
fprintf(fid, '.best { background-color: yellow; font-weight: bold; }\n');
fprintf(fid, '.note { font-size: 10pt; text-align: center; }\n');
fprintf(fid, '</style>\n');

fprintf(fid, '</head>\n<body>\n');

fprintf(fid, '<h1>All Methods Summary</h1>\n');
fprintf(fid, '<p><b>Train ratio:</b> %.0f%% &nbsp;&nbsp; <b>Monte Carlo runs:</b> %d</p>\n', train_ratio * 100, n_mc);
fprintf(fid, '<p>Rows are method families and columns are datasets. For each metric and aggregation rule, the best value in each dataset column is highlighted. Lower values are better for SMSE, RMSE, train time, and test time.</p>\n');

for met = 1:numel(metrics_names)
    metric_name = metrics_names{met};

    fprintf(fid, '\n<h1>Metric: %s &nbsp; (Train=%.0f%%, MC=%d)</h1>\n', ...
        html_escape(metric_name), train_ratio * 100, n_mc);

    for ai = 1:numel(agg_list)
        agg = agg_list{ai};

        [cell_text, best_mask] = build_table_values( ...
            mean_results, std_results, datasets, methods_names, ...
            met, agg, row_families);

        fprintf(fid, '\n<h2>%s</h2>\n', html_escape(agg));
        fprintf(fid, '<table>\n');

        % First header row: aggregation name spanning dataset columns.
        fprintf(fid, '<tr>\n');
        fprintf(fid, '<th style="width:18%%;"></th>\n');
        fprintf(fid, '<th class="agg" colspan="%d">%s</th>\n', numel(dataset_labels), html_escape(agg));
        fprintf(fid, '</tr>\n');

        % Second header row: dataset labels.
        fprintf(fid, '<tr>\n');
        fprintf(fid, '<th>Method family</th>\n');
        for d = 1:numel(dataset_labels)
            fprintf(fid, '<th>%s</th>\n', html_escape(dataset_labels{d}));
        end
        fprintf(fid, '</tr>\n');

        % Data rows.
        for r = 1:numel(row_families)
            fprintf(fid, '<tr>\n');
            fprintf(fid, '<td class="rowhead">%s</td>\n', html_escape(row_families{r}));

            for d = 1:numel(dataset_labels)
                txt = html_escape(cell_text{r,d});

                if best_mask(r,d) && ~strcmp(cell_text{r,d}, '-')
                    fprintf(fid, '<td class="best">%s</td>\n', txt);
                else
                    fprintf(fid, '<td>%s</td>\n', txt);
                end
            end

            fprintf(fid, '</tr>\n');
        end

        fprintf(fid, '</table>\n');
    end
end

fprintf(fid, '</body>\n</html>\n');
fclose(fid);

fprintf('\nEditable HTML tables saved to:\n%s\n', HtmlPath);
fprintf('Open this file in a browser or Microsoft Word, then copy/paste the editable tables into your report.\n');

%% ========================================================================
%% Build table values
%% ========================================================================

function [cell_text, best_mask] = build_table_values( ...
    mean_results, std_results, datasets, methods_names, ...
    met, agg, row_families)

    n_rows = numel(row_families);
    n_cols = numel(datasets);

    value_mean = NaN(n_rows, n_cols);
    cell_text  = cell(n_rows, n_cols);

    for r = 1:n_rows
        family = row_families{r};
        target_name = sprintf('%s-%s', family, agg);
        mi = find(strcmp(methods_names, target_name), 1);

        for d = 1:n_cols
            if isempty(mi)
                cell_text{r,d} = '-';
            else
                mv = mean_results(d, mi, met);
                sv = std_results(d,  mi, met);

                value_mean(r,d) = mv;
                cell_text{r,d}  = format_result_value(mv, sv, met);
            end
        end
    end

    best_mask = false(n_rows, n_cols);

    for d = 1:n_cols
        col_vals = value_mean(:, d);
        valid = isfinite(col_vals);

        % For communication / iteration metrics, ignore zero because it is
        % displayed as "-".
        if met >= 5
            valid = valid & col_vals > 0;
        end

        if any(valid)
            best_val = min(col_vals(valid));
            best_mask(:, d) = valid & abs(col_vals - best_val) < 1e-12;
        end
    end
end

%% ========================================================================
%% Formatting helpers
%% ========================================================================

function value_text = format_result_value(mv, sv, met)
    if isnan(mv)
        value_text = '-';
        return;
    end

    if met >= 5
        % Communication and iteration columns: integer display.
        % 0 is shown as '-'.
        if mv == 0 || isnan(mv)
            value_text = '-';
        else
            value_text = sprintf('%.0f', mv);
        end
    elseif met == 3 || met == 4
        value_text = sprintf('%.2f±%.2f', mv, sv);
    else
        value_text = sprintf('%.4f±%.4f', mv, sv);
    end
end

function out = html_escape(in)
    if isnumeric(in)
        in = num2str(in);
    end
    out = char(in);
    out = strrep(out, '&', '&amp;');
    out = strrep(out, '<', '&lt;');
    out = strrep(out, '>', '&gt;');
end
