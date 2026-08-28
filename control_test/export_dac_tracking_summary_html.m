function export_dac_tracking_summary_html(method_list, form_tag, result_folder, output_folder)
%EXPORT_DAC_TRACKING_SUMMARY_HTML Export DAC control summary as Word-friendly HTML.
%
% Usage:
%   export_dac_tracking_summary_html
%   export_dac_tracking_summary_html({'poe','gpoe','moe','bcm','rbcm'}, 'formation')

if nargin < 1 || isempty(method_list)
    method_list = {'poe','gpoe','moe','bcm','rbcm'};
end
if ischar(method_list) || isstring(method_list)
    method_list = cellstr(method_list);
end
if nargin < 2 || isempty(form_tag)
    form_tag = 'formation';
end
if nargin < 3 || isempty(result_folder)
    result_folder = fullfile('result', 'Diagnostics');
end
if nargin < 4 || isempty(output_folder)
    output_folder = fullfile('result', 'Figures');
end
if ~isfolder(output_folder)
    mkdir(output_folder);
end

rows = struct('method', {}, 'final_error', {}, 'trigger_avg', {}, ...
    'trigger_per_agent', {}, 'file_name', {});
for method_idx = 1:numel(method_list)
    method = lower(method_list{method_idx});
    result_file = find_dac_result_file(result_folder, method, form_tag);
    if isempty(result_file)
        warning('No DAC result file found for method "%s".', method);
        continue;
    end

    data = load(result_file, 'TrackingError_vector', 'dac_broadcast_count_set');
    trigger_per_agent = sum(data.dac_broadcast_count_set, 2);
    rows(end+1).method = upper(method); %#ok<AGROW>
    rows(end).final_error = data.TrackingError_vector(end);
    rows(end).trigger_avg = mean(trigger_per_agent);
    rows(end).trigger_per_agent = trigger_per_agent(:).';
    [~, name, ext] = fileparts(result_file);
    rows(end).file_name = [name, ext];
end

html_file = fullfile(output_folder, sprintf('DAC_Tracking_Summary_%s.html', form_tag));
fid = fopen(html_file, 'w');
if fid < 0
    error('Cannot open output file: %s', html_file);
end
cleanup_obj = onCleanup(@() fclose(fid));

fprintf(fid, '<!DOCTYPE html>\n<html><head><meta charset="UTF-8">\n');
fprintf(fid, '<style>\n');
fprintf(fid, 'body{font-family:"Times New Roman",serif;margin:24px;color:#111;}\n');
fprintf(fid, 'table{border-collapse:collapse;font-size:12pt;}\n');
fprintf(fid, 'th,td{border:1px solid #555;padding:6px 10px;text-align:center;}\n');
fprintf(fid, 'th{background:#e9eef7;font-weight:bold;}\n');
fprintf(fid, 'td.method{font-weight:bold;background:#f7f7f7;}\n');
fprintf(fid, 'caption{caption-side:top;font-weight:bold;font-size:14pt;margin-bottom:8px;}\n');
fprintf(fid, '</style></head><body>\n');
fprintf(fid, '<table>\n');
fprintf(fid, '<caption>IP-DAC Final Tracking Error and Broadcast Triggers [%s]</caption>\n', form_tag);
fprintf(fid, '<tr><th>Method</th><th>Final tracking error</th><th>Avg. DAC broadcasts / agent</th><th>Broadcasts per agent</th></tr>\n');
for idx = 1:numel(rows)
    fprintf(fid, '<tr>');
    fprintf(fid, '<td class="method">%s</td>', rows(idx).method);
    fprintf(fid, '<td>%.4f</td>', rows(idx).final_error);
    fprintf(fid, '<td>%.2f</td>', rows(idx).trigger_avg);
    fprintf(fid, '<td>%s</td>', format_vector(rows(idx).trigger_per_agent));
    fprintf(fid, '</tr>\n');
end
fprintf(fid, '</table>\n');
fprintf(fid, '<p style="font-size:10pt;margin-top:12px;">Trigger count is reported at the agent broadcast level over physical simulation time. ');
fprintf(fid, 'DAC uses Kia et al. law (17) with epsilon_i/(2 sqrt(d_i^out)) = 0.04.</p>\n');
fprintf(fid, '</body></html>\n');

fprintf('Wrote DAC summary HTML: %s\n', html_file);
end

function result_file = find_dac_result_file(result_folder, method, form_tag)
preferred = fullfile(result_folder, ...
    sprintf('diag_%s_law17_fig3eps004_unknown0p20_dist0p10_%s.mat', method, form_tag));
if isfile(preferred)
    result_file = preferred;
    return;
end

files = dir(fullfile(result_folder, sprintf('diag_%s*%s*.mat', method, form_tag)));
files = files(~contains({files.name}, '_ac_'));
if isempty(files)
    result_file = '';
    return;
end

[~, order] = sort([files.datenum], 'descend');
files = files(order);
for idx = 1:numel(files)
    candidate = fullfile(files(idx).folder, files(idx).name);
    names = string({whos('-file', candidate).name});
    if any(names == "TrackingError_vector") && any(names == "dac_broadcast_count_set")
        result_file = candidate;
        return;
    end
end
result_file = '';
end

function text = format_vector(values)
parts = arrayfun(@(x) sprintf('%.0f', x), values, 'UniformOutput', false);
text = strjoin(parts, ', ');
end
