%% export_summary_tables_to_html_dataset_columns.m
% One editable HTML table per metric.
%
% Layout:
%   Rows:
%       Aggregation | Method family
%   Columns:
%       KIN40K | POL | PUMA | SARCOS
%
% This layout is Word-friendly and keeps one metric as one large table.

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

agg_list = {'MOE','POE','GPOE','BCM','RBCM'};
row_families = {'LoG','CEN','IP-DAC','IP-AC','TP-DAC','TP-AC','NBR'};

dataset_labels = datasets;
dataset_labels(strcmp(dataset_labels, 'PUMADYN32NM')) = {'PUMA'};

OutFolder = fullfile(ProjectRoot, 'Result', 'Dataset');
if ~exist(OutFolder, 'dir')
    mkdir(OutFolder);
end

HtmlPath = fullfile(OutFolder, 'All_Methods_Summary_dataset_columns.html');

fid = fopen(HtmlPath, 'w');
if fid < 0
    error('Cannot open HTML file for writing:\n%s', HtmlPath);
end

fprintf(fid, '<!DOCTYPE html>\n');
fprintf(fid, '<html>\n<head>\n<meta charset="UTF-8">\n');
fprintf(fid, '<title>All Methods Summary</title>\n');

fprintf(fid, '<style>\n');
fprintf(fid, 'body { font-family: "Times New Roman", serif; margin: 20px; }\n');
fprintf(fid, 'h1 { color: #b00000; font-size: 20pt; text-align: center; margin-top: 24px; }\n');
fprintf(fid, 'p { font-size: 10.5pt; }\n');
fprintf(fid, 'table { border-collapse: collapse; margin: 12px auto 30px auto; width: 92%%; table-layout: fixed; }\n');
fprintf(fid, 'th, td { border: 1px solid black; padding: 5px 7px; text-align: center; font-size: 10pt; }\n');
fprintf(fid, 'th { background-color: #cfe8f3; font-weight: bold; }\n');
fprintf(fid, '.aggcol { font-weight: bold; background-color: #f2f2f2; }\n');
fprintf(fid, '.methodcol { font-weight: bold; background-color: #f7f7f7; }\n');
fprintf(fid, '.ipdac { font-weight: bold; background-color: #fff7cc; }\n');
fprintf(fid, '.missing { color: #777777; }\n');
fprintf(fid, '@page { size: landscape; margin: 0.6in; }\n');
fprintf(fid, '</style>\n');

fprintf(fid, '</head>\n<body>\n');

fprintf(fid, '<h1>All Methods Summary</h1>\n');
fprintf(fid, '<p><b>Train ratio:</b> %.0f%% &nbsp;&nbsp; <b>Monte Carlo runs:</b> %d</p>\n', ...
    train_ratio * 100, n_mc);

fprintf(fid, '<p>Each metric is shown as one table. Rows are grouped by aggregation rule, and columns are datasets. The IP-DAC row is lightly emphasized for comparison.</p>\n');

for met = 1:numel(metrics_names)
    metric_name = metrics_names{met};

    fprintf(fid, '\n<h1>Metric: %s &nbsp; (Train=%.0f%%, MC=%d)</h1>\n', ...
        html_escape(metric_name), train_ratio * 100, n_mc);

    fprintf(fid, '<table>\n');

    fprintf(fid, '<tr>\n');
    fprintf(fid, '<th style="width:12%%;">Aggregation</th>\n');
    fprintf(fid, '<th style="width:12%%;">Method</th>\n');
    for d = 1:numel(dataset_labels)
        fprintf(fid, '<th>%s</th>\n', html_escape(dataset_labels{d}));
    end
    fprintf(fid, '</tr>\n');

    for a = 1:numel(agg_list)
        agg = agg_list{a};

        for r = 1:numel(row_families)
            family = row_families{r};
            target_name = sprintf('%s-%s', family, agg);
            mi = find(strcmp(methods_names, target_name), 1);

            if strcmp(family, 'IP-DAC')
                row_class = ' class="ipdac"';
            else
                row_class = '';
            end

            fprintf(fid, '<tr%s>\n', row_class);

            if r == 1
                fprintf(fid, '<td class="aggcol" rowspan="%d">%s</td>\n', ...
                    numel(row_families), html_escape(agg));
            end

            fprintf(fid, '<td class="methodcol">%s</td>\n', html_escape(family));

            for d = 1:numel(datasets)
                if isempty(mi)
                    txt_raw = '-';
                else
                    mv = mean_results(d, mi, met);
                    sv = std_results(d,  mi, met);
                    txt_raw = format_result_value(mv, sv, met);
                end

                txt = html_escape(txt_raw);

                if strcmp(txt_raw, '-')
                    fprintf(fid, '<td class="missing">-</td>\n');
                else
                    fprintf(fid, '<td>%s</td>\n', txt);
                end
            end

            fprintf(fid, '</tr>\n');
        end
    end

    fprintf(fid, '</table>\n');
end

fprintf(fid, '</body>\n</html>\n');
fclose(fid);

fprintf('\nWord-friendly dataset-column HTML saved to:\n%s\n', HtmlPath);
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