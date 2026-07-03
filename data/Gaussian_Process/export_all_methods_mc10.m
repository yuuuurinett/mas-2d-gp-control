%% export_dataset_columns_with_LMH_tables.m
% Export dataset-column summary tables with Low/Medium/High rows.
%
% Layout:
%   One table per metric.
%
% Columns:
%   Aggregation | Method | M Level | KIN40K | POL | PUMA | SARCOS
%
% Rows:
%   Aggregation blocks: MOE / POE / GPOE / BCM / RBCM
%   Method rows:
%       LoG
%       CEN
%       IP-DAC Low / Medium / High
%       IP-AC  Low / Medium / High
%       TP-DAC
%       TP-AC
%       NBR
%
% For normal metrics:
%   IP-DAC/IP-AC read Low/Medium/High from M-ablation files.
%
% For Trigger_Ratio_Train(%):
%   IP-DAC/IP-AC show only one final-MC row.
%   Other methods show "-".
%
% Outputs:
%   Result/Dataset/All_Methods_Summary_dataset_columns_with_LMH.html
%   Result/Dataset/All_Methods_Summary_dataset_columns_with_LMH.md

clear; clc;

ProjectRoot = fileparts(mfilename('fullpath'));
cd(ProjectRoot);
addpath(genpath(ProjectRoot));

DatasetResultRoot = fullfile(ProjectRoot, 'Result', 'Dataset');
SummaryPath = fullfile(DatasetResultRoot, 'All_Methods_Summary.mat');

if ~isfile(SummaryPath)
    error('Cannot find summary file:\n%s', SummaryPath);
end

load(SummaryPath, ...
    'mean_results', 'std_results', ...
    'datasets', 'methods_names', 'metrics_names', ...
    'train_ratio', 'n_mc', 'tr_tag');

%% ===================== Config =====================

% Keep the same order as your screenshot if you prefer:
agg_list = {'MOE','POE','GPOE','BCM','RBCM'};

agg_to_file = containers.Map( ...
    {'MOE','GPOE','POE','BCM','RBCM'}, ...
    {'moe','gpoe','poe','bcm','rbcm'});

method_order = {'LoG','CEN','IP-DAC','IP-AC','TP-DAC','TP-AC','NBR'};
ip_methods = {'IP-DAC','IP-AC'};

% Low / Medium / High for M-ablation rows
M_levels_default = [500, 1500, 2500];   % KIN40K / POL / SARCOS
M_levels_puma    = [100, 300, 500];     % PUMADYN32NM
M_level_names    = {'Low','Medium','High'};

HtmlPath = fullfile(DatasetResultRoot, 'All_Methods_Summary_dataset_columns_with_LMH.html');
MdPath   = fullfile(DatasetResultRoot, 'All_Methods_Summary_dataset_columns_with_LMH.md');

%% ===================== Write HTML =====================

fid = fopen(HtmlPath, 'w');
if fid < 0
    error('Cannot open HTML file:\n%s', HtmlPath);
end

fprintf(fid, '<!DOCTYPE html>\n<html>\n<head>\n<meta charset="UTF-8">\n');
fprintf(fid, '<title>Dataset Results Summary with Low/Medium/High</title>\n');

fprintf(fid, '<style>\n');
fprintf(fid, '@page { size: A4 landscape; margin: 0.45in; }\n');
fprintf(fid, 'body { font-family: "Times New Roman", serif; margin: 12px; }\n');
fprintf(fid, 'h1 { color: #2f6fb3; font-size: 18pt; margin: 8px 0 12px 0; }\n');
fprintf(fid, 'h2 { color: #2f6fb3; font-size: 17pt; margin: 16px 0 8px 0; }\n');
fprintf(fid, 'p { font-size: 9pt; margin: 3px 0; }\n');
fprintf(fid, 'table { border-collapse: collapse; width: 96%%; table-layout: fixed; margin: 8px 0 22px 0; page-break-inside: avoid; }\n');
fprintf(fid, 'th, td { border: 1px solid black; padding: 3px 5px; text-align: center; vertical-align: middle; font-size: 8.5pt; line-height: 1.15; }\n');
fprintf(fid, 'th { background-color: #cfe8f3; font-weight: bold; }\n');
fprintf(fid, '.agg { font-weight: bold; background-color: #f2f2f2; width: 9%%; }\n');
fprintf(fid, '.method { font-weight: bold; background-color: #f7f7f7; width: 8%%; }\n');
fprintf(fid, '.level { font-weight: bold; background-color: #fbfbfb; width: 7%%; }\n');
fprintf(fid, '.ip { font-weight: bold; }\n');
fprintf(fid, '.missing { color: #777777; }\n');
fprintf(fid, '.metric-block { page-break-after: always; }\n');
fprintf(fid, '</style>\n');

fprintf(fid, '</head>\n<body>\n');

fprintf(fid, '<h1>Dataset Results Summary</h1>\n');
fprintf(fid, '<p><b>Train ratio:</b> %.0f%% &nbsp;&nbsp; <b>MC:</b> %d</p>\n', train_ratio * 100, n_mc);
fprintf(fid, '<p>For IP-DAC/IP-AC, Low/Medium/High rows are shown for accuracy/time metrics.</p>\n');
fprintf(fid, '<p>KIN40K/POL/SARCOS: Low=500, Medium=1500, High=2500. PUMADYN32NM: Low=100, Medium=300, High=500.</p>\n');
fprintf(fid, '<p>For Trigger_Ratio_Train(%%), IP-DAC/IP-AC show one final-MC value only.</p>\n');

for met = 1:numel(metrics_names)
    metric_name = metrics_names{met};

    fprintf(fid, '<div class="metric-block">\n');
    fprintf(fid, '<h2>%s__(Train=%.0f%%, MC=%d)</h2>\n', ...
        escape_html(metric_name), train_ratio * 100, n_mc);

    fprintf(fid, '<table>\n');

    % Header
    fprintf(fid, '<tr>\n');
    fprintf(fid, '<th>Aggregation</th>');
    fprintf(fid, '<th>Method</th>');
    fprintf(fid, '<th>M Level</th>');
    for d = 1:numel(datasets)
        fprintf(fid, '<th>%s</th>', escape_html(short_dataset_name(datasets{d})));
    end
    fprintf(fid, '\n</tr>\n');

    for a_idx = 1:numel(agg_list)
        agg = agg_list{a_idx};

        % Count rows for rowspan.
        % Normal metrics: LoG,CEN,IP-DACx3,IP-ACx3,TP-DAC,TP-AC,NBR = 11 rows
        % Trigger ratio: all methods single row = 7 rows
        if met == 7
            agg_rowspan = 7;
        else
            agg_rowspan = 11;
        end

        agg_cell_printed = false;

        for mp = 1:numel(method_order)
            method_prefix = method_order{mp};
            is_ip = ismember(method_prefix, ip_methods);

            if is_ip && met ~= 7
                % IP-DAC/IP-AC: Low / Medium / High rows
                for lv = 1:numel(M_level_names)
                    level_name = M_level_names{lv};

                    fprintf(fid, '<tr>\n');

                    if ~agg_cell_printed
                        fprintf(fid, '<td rowspan="%d" class="agg">%s</td>', ...
                            agg_rowspan, escape_html(agg));
                        agg_cell_printed = true;
                    end

                    if lv == 1
                        fprintf(fid, '<td rowspan="3" class="method ip">%s</td>', escape_html(method_prefix));
                    end

                    fprintf(fid, '<td class="level">%s</td>', escape_html(level_name));

                    for d = 1:numel(datasets)
                        dname = datasets{d};

                        if strcmpi(dname, 'PUMADYN32NM')
                            M_val = M_levels_puma(lv);
                        else
                            M_val = M_levels_default(lv);
                        end

                        agg_lower = agg_to_file(agg);
                        vals = read_ip_m_ablation_values( ...
                            DatasetResultRoot, dname, agg_lower, M_val, tr_tag, n_mc, met);

                        mv = mean(vals, 'omitnan');
                        sv = std(vals, 0, 'omitnan');
                        value_text = format_table_value(mv, sv, met);

                        fprintf(fid, '<td>%s</td>', escape_html(value_text));
                    end

                    fprintf(fid, '\n</tr>\n');
                end

            else
                % Single row for non-IP methods, and also for Trigger_Ratio_Train
                fprintf(fid, '<tr>\n');

                if ~agg_cell_printed
                    fprintf(fid, '<td rowspan="%d" class="agg">%s</td>', ...
                        agg_rowspan, escape_html(agg));
                    agg_cell_printed = true;
                end

                fprintf(fid, '<td class="method">%s</td>', escape_html(method_prefix));
                fprintf(fid, '<td class="level">-</td>');

                for d = 1:numel(datasets)
                    dname = datasets{d};
                    target_name = sprintf('%s-%s', method_prefix, agg);
                    mi = find(strcmp(methods_names, target_name), 1);

                    if isempty(mi)
                        value_text = '-';
                    else
                        if (met == 5 || met == 7) && ~is_ip
                            value_text = '-';
                        else
                            mv = mean_results(d, mi, met);
                            sv = std_results(d, mi, met);
                            value_text = format_table_value(mv, sv, met);
                        end
                    end

                    fprintf(fid, '<td>%s</td>', escape_html(value_text));
                end

                fprintf(fid, '\n</tr>\n');
            end
        end
    end

    fprintf(fid, '</table>\n');
    fprintf(fid, '</div>\n');
end

fprintf(fid, '</body>\n</html>\n');
fclose(fid);

fprintf('HTML saved:\n%s\n', HtmlPath);

%% ===================== Write Markdown =====================

fid = fopen(MdPath, 'w');
if fid < 0
    error('Cannot open Markdown file:\n%s', MdPath);
end

fprintf(fid, '# Dataset Results Summary with Low/Medium/High\n\n');
fprintf(fid, '**Train ratio:** %.0f%%  \n', train_ratio * 100);
fprintf(fid, '**MC:** %d  \n\n', n_mc);
fprintf(fid, 'KIN40K/POL/SARCOS: Low=500, Medium=1500, High=2500.  \n');
fprintf(fid, 'PUMADYN32NM: Low=100, Medium=300, High=500.  \n\n');

for met = 1:numel(metrics_names)
    metric_name = metrics_names{met};

    fprintf(fid, '\n\n# %s__(Train=%.0f%%, MC=%d)\n\n', ...
        metric_name, train_ratio * 100, n_mc);

    fprintf(fid, '| Aggregation | Method | M Level |');
    for d = 1:numel(datasets)
        fprintf(fid, ' %s |', short_dataset_name(datasets{d}));
    end
    fprintf(fid, '\n');

    fprintf(fid, '|---|---|---|');
    for d = 1:numel(datasets)
        fprintf(fid, '---:|');
    end
    fprintf(fid, '\n');

    for a_idx = 1:numel(agg_list)
        agg = agg_list{a_idx};

        for mp = 1:numel(method_order)
            method_prefix = method_order{mp};
            is_ip = ismember(method_prefix, ip_methods);

            if is_ip && met ~= 7
                for lv = 1:numel(M_level_names)
                    level_name = M_level_names{lv};

                    if lv == 1
                        agg_cell = agg;
                        method_cell = method_prefix;
                    else
                        agg_cell = '';
                        method_cell = '';
                    end

                    fprintf(fid, '| %s | %s | %s |', agg_cell, method_cell, level_name);

                    for d = 1:numel(datasets)
                        dname = datasets{d};

                        if strcmpi(dname, 'PUMADYN32NM')
                            M_val = M_levels_puma(lv);
                        else
                            M_val = M_levels_default(lv);
                        end

                        agg_lower = agg_to_file(agg);
                        vals = read_ip_m_ablation_values( ...
                            DatasetResultRoot, dname, agg_lower, M_val, tr_tag, n_mc, met);

                        mv = mean(vals, 'omitnan');
                        sv = std(vals, 0, 'omitnan');
                        value_text = format_table_value(mv, sv, met);

                        fprintf(fid, ' %s |', value_text);
                    end

                    fprintf(fid, '\n');
                end

            else
                fprintf(fid, '| %s | %s | - |', agg, method_prefix);

                for d = 1:numel(datasets)
                    dname = datasets{d};
                    target_name = sprintf('%s-%s', method_prefix, agg);
                    mi = find(strcmp(methods_names, target_name), 1);

                    if isempty(mi)
                        value_text = '-';
                    else
                        if (met == 5 || met == 7) && ~is_ip
                            value_text = '-';
                        else
                            mv = mean_results(d, mi, met);
                            sv = std_results(d, mi, met);
                            value_text = format_table_value(mv, sv, met);
                        end
                    end

                    fprintf(fid, ' %s |', value_text);
                end

                fprintf(fid, '\n');
            end
        end
    end
end

fclose(fid);

fprintf('Markdown saved:\n%s\n', MdPath);

%% ========================================================================
%% Helper functions
%% ========================================================================

function vals = read_ip_m_ablation_values(DatasetResultRoot, dname, agg_lower, M_val, tr_tag, n_mc, met)
    vals = NaN(1, n_mc);

    for mc = 1:n_mc
        fname = sprintf('%s_M%d_tr%d_mc%d.mat', agg_lower, M_val, tr_tag, mc);
        fpath = fullfile(DatasetResultRoot, dname, fname);

        if ~isfile(fpath)
            continue;
        end

        try
            res = load(fpath);

            switch met
                case 1
                    if isfield(res,'smse')
                        vals(mc) = res.smse;
                    end

                case 2
                    if isfield(res,'rmse')
                        vals(mc) = res.rmse;
                    end

                case 3
                    if isfield(res,'t_train_per_point')
                        vals(mc) = res.t_train_per_point;
                    end

                case 4
                    if isfield(res,'t_test_per_point')
                        vals(mc) = res.t_test_per_point;
                    end

                case 5
                    if isfield(res,'comm_train')
                        vals(mc) = res.comm_train;
                    elseif isfield(res,'Trigger_Times')
                        vals(mc) = res.Trigger_Times;
                    end

                case 6
                    if isfield(res,'comm_test')
                        vals(mc) = res.comm_test;
                    end
            end
        catch
            vals(mc) = NaN;
        end
    end
end

function value_text = format_table_value(mv, sv, met)
    if isnan(mv)
        value_text = '-';
        return;
    end

    if met == 7
        value_text = sprintf('%.1f±%.1f%%', mv, sv);

    elseif met == 6
        if mv == 0
            value_text = '-';
        else
            value_text = sprintf('%.0f', mv);
        end

    elseif met == 5
        value_text = sprintf('%.0f', mv);

    elseif met == 3 || met == 4
        value_text = sprintf('%.3f±%.3f', mv, sv);

    else
        value_text = sprintf('%.4f±%.4f', mv, sv);
    end
end

function name = short_dataset_name(name)
    if strcmpi(name, 'PUMADYN32NM')
        name = 'PUMA';
    end
end

function s = escape_html(s)
    s = char(s);
    s = strrep(s, '&', '&amp;');
    s = strrep(s, '<', '&lt;');
    s = strrep(s, '>', '&gt;');
    s = strrep(s, '"', '&quot;');
end