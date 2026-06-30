%% export_summary_tables_to_html_big_metric_tables.m
% One editable HTML table per metric, without automatic highlighting.
%
% Table layout:
%   Columns: Aggregation groups = MOE / POE / GPOE / BCM / RBCM
%            Each aggregation group contains datasets:
%            KIN40K / POL / PUMA / SARCOS
%
%   Rows:    LoG / CEN / IP-DAC / IP-AC / TP-DAC / TP-AC / NBR
%
% No automatic best-value highlighting.
% The purpose is to compare IP-DAC against other method families manually.

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

% Horizontal aggregation group order
agg_list = {'MOE','POE','GPOE','BCM','RBCM'};

% Vertical method-family order
row_families = {'LoG','CEN','IP-DAC','IP-AC','TP-DAC','TP-AC','NBR'};

% Compact dataset labels
dataset_labels = datasets;
dataset_labels(strcmp(dataset_labels, 'PUMADYN32NM')) = {'PUMA'};

OutFolder = fullfile(ProjectRoot, 'Result', 'Dataset');
if ~exist(OutFolder, 'dir')
    mkdir(OutFolder);
end

HtmlPath = fullfile(OutFolder, 'All_Methods_Summary_big_metric_tables.html');

fid = fopen(HtmlPath, 'w');
if fid < 0
    error('Cannot open HTML file for writing:\n%s', HtmlPath);
end

%% ===================== HTML header =====================

fprintf(fid, '<!DOCTYPE html>\n');
fprintf(fid, '<html>\n<head>\n<meta charset="UTF-8">\n');
fprintf(fid, '<title>All Methods Summary</title>\n');

fprintf(fid, '<style>\n');
fprintf(fid, 'body { font-family: "Times New Roman", serif; margin: 20px; }\n');
fprintf(fid, 'h1 { color: #b00000; font-size: 20pt; text-align: center; margin-top: 26px; }\n');
fprintf(fid, 'p { font-size: 10.5pt; }\n');

% Wide table, so use compact font size.
fprintf(fid, 'table { border-collapse: collapse; margin: 12px auto 30px auto; width: 100%%; table-layout: fixed; }\n');
fprintf(fid, 'th, td { border: 1px solid black; padding: 4px 4px; text-align: center; font-size: 8.5pt; }\n');
fprintf(fid, 'th { background-color: #cfe8f3; font-weight: bold; }\n');
fprintf(fid, '.rowhead { font-weight: bold; background-color: #f7f7f7; width: 7%%; }\n');
fprintf(fid, '.agghead { background-color: #cfe8f3; font-weight: bold; font-size: 10pt; }\n');
fprintf(fid, '.datasethead { background-color: #e8f4fa; font-weight: bold; }\n');
fprintf(fid, '.missing { color: #777777; }\n');

% Word landscape layout hint
fprintf(fid, '@page { size: landscape; margin: 0.5in; }\n');

fprintf(fid, '</style>\n');
fprintf(fid, '</head>\n<body>\n');

fprintf(fid, '<h1>All Methods Summary</h1>\n');
fprintf(fid, '<p><b>Train ratio:</b> %.0f%% &nbsp;&nbsp; <b>Monte Carlo runs:</b> %d</p>\n', ...
    train_ratio * 100, n_mc);

fprintf(fid, '<p>Each metric is shown as one large table. Columns are grouped by aggregation rule, and each aggregation group contains the four datasets. Rows are method families. No automatic best-value highlighting is used, so IP-DAC can be compared manually against the other method families.</p>\n');

%% ===================== One big table per metric =====================

for met = 1:numel(metrics_names)
    metric_name = metrics_names{met};

    fprintf(fid, '\n<h1>Metric: %s &nbsp; (Train=%.0f%%, MC=%d)</h1>\n', ...
        html_escape(metric_name), train_ratio * 100, n_mc);

    n_rows = numel(row_families);
    n_agg  = numel(agg_list);
    n_data = numel(datasets);

    value_text = cell(n_rows, n_agg, n_data);

    for r = 1:n_rows
        family = row_families{r};

        for a = 1:n_agg
            agg = agg_list{a};
            target_name = sprintf('%s-%s', family, agg);
            mi = find(strcmp(methods_names, target_name), 1);

            for d = 1:n_data
                if isempty(mi)
                    value_text{r,a,d} = '-';
                else
                    mv = mean_results(d, mi, met);
                    sv = std_results(d,  mi, met);
                    value_text{r,a,d} = format_result_value(mv, sv, met);
                end
            end
        end
    end

    %% --------------------- Write HTML table ---------------------

    fprintf(fid, '<table>\n');

    % First header row: aggregation names
    fprintf(fid, '<tr>\n');
    fprintf(fid, '<th class="rowhead"></th>\n');

    for a = 1:n_agg
        fprintf(fid, '<th class="agghead" colspan="%d">%s</th>\n', ...
            n_data, html_escape(agg_list{a}));
    end

    fprintf(fid, '</tr>\n');

    % Second header row: dataset names under each aggregation
    fprintf(fid, '<tr>\n');
    fprintf(fid, '<th class="rowhead">Method</th>\n');

    for a = 1:n_agg
        for d = 1:n_data
            fprintf(fid, '<th class="datasethead">%s</th>\n', ...
                html_escape(dataset_labels{d}));
        end
    end

    fprintf(fid, '</tr>\n');

    % Data rows
    for r = 1:n_rows
        fprintf(fid, '<tr>\n');
        fprintf(fid, '<td class="rowhead">%s</td>\n', html_escape(row_families{r}));

        for a = 1:n_agg
            for d = 1:n_data
                txt_raw = value_text{r,a,d};
                txt = html_escape(txt_raw);

                if strcmp(txt_raw, '-')
                    fprintf(fid, '<td class="missing">-</td>\n');
                else
                    fprintf(fid, '<td>%s</td>\n', txt);
                end
            end
        end

        fprintf(fid, '</tr>\n');
    end

    fprintf(fid, '</table>\n');
end

fprintf(fid, '</body>\n</html>\n');
fclose(fid);

fprintf('\nBig editable HTML tables saved to:\n%s\n', HtmlPath);
fprintf('Open this HTML file with browser or Microsoft Word, then copy/paste into your report.\n');

%% ========================================================================
%% Helper functions
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
        % Time
        value_text = sprintf('%.2f±%.2f', mv, sv);
    else
        % SMSE / RMSE / MSLL etc.
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