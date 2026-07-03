%% generate_summary_html.m
% Generates a clean HTML summary table from M-ablation .mat files.
% Output: Result/Dataset/Summary_Table.html
% Open in browser, then copy-paste into Word.

clear; clc;

ProjectRoot = fileparts(mfilename('fullpath'));
cd(ProjectRoot);

%% ===================== Configuration =====================

datasets     = {'KIN40K','POL','PUMADYN32NM','SARCOS'};
agg_list     = {'MOE','GPOE','POE','BCM','RBCM'};
col_prefixes = {'LoG','CEN','IP-DAC','IP-AC','TP-DAC','TP-AC','NBR'};

train_ratio = 0.4;
tr_tag      = round(train_ratio * 100);
n_mc        = 10;
mc_list     = 1:n_mc;

% M levels
M_default    = 100:100:2500;
M_puma       = 100:20:500;
LMH_default  = [500, 1500, 2500];
LMH_puma     = [100, 300,  500];

% File patterns (same as run_mc_dataset.m)
methods_dict = {
    'LoG-MOE',      'log_moe_mc%d.mat',        false;
    'LoG-GPOE',     'log_gpoe_mc%d.mat',        false;
    'LoG-POE',      'log_poe_mc%d.mat',         false;
    'LoG-BCM',      'log_bcm_mc%d.mat',         false;
    'LoG-RBCM',     'log_rbcm_mc%d.mat',        false;
    'IP-DAC-MOE',   'moe_tr%d_mc%d.mat',        true;
    'IP-DAC-GPOE',  'gpoe_tr%d_mc%d.mat',       true;
    'IP-DAC-POE',   'poe_tr%d_mc%d.mat',        true;
    'IP-DAC-BCM',   'bcm_tr%d_mc%d.mat',        true;
    'IP-DAC-RBCM',  'rbcm_tr%d_mc%d.mat',       true;
    'IP-AC-MOE',    'moe_ac_tr%d_mc%d.mat',     true;
    'IP-AC-GPOE',   'gpoe_ac_tr%d_mc%d.mat',    true;
    'IP-AC-POE',    'poe_ac_tr%d_mc%d.mat',     true;
    'IP-AC-BCM',    'bcm_ac_tr%d_mc%d.mat',     true;
    'IP-AC-RBCM',   'rbcm_ac_tr%d_mc%d.mat',   true;
    'TP-DAC-MOE',   'moe_tp_tr%d_mc%d.mat',     true;
    'TP-DAC-GPOE',  'gpoe_tp_tr%d_mc%d.mat',    true;
    'TP-DAC-POE',   'poe_tp_tr%d_mc%d.mat',     true;
    'TP-DAC-BCM',   'bcm_tp_tr%d_mc%d.mat',     true;
    'TP-DAC-RBCM',  'rbcm_tp_tr%d_mc%d.mat',   true;
    'CEN-MOE',      'moe_cen_tr%d_mc%d.mat',    true;
    'CEN-GPOE',     'gpoe_cen_tr%d_mc%d.mat',   true;
    'CEN-POE',      'poe_cen_tr%d_mc%d.mat',    true;
    'CEN-BCM',      'bcm_cen_tr%d_mc%d.mat',    true;
    'CEN-RBCM',     'rbcm_cen_tr%d_mc%d.mat',  true;
    'NBR-MOE',      'moe_nbr_tr%d_mc%d.mat',    true;
    'NBR-GPOE',     'gpoe_nbr_tr%d_mc%d.mat',   true;
    'NBR-POE',      'poe_nbr_tr%d_mc%d.mat',    true;
    'NBR-BCM',      'bcm_nbr_tr%d_mc%d.mat',    true;
    'NBR-RBCM',     'rbcm_nbr_tr%d_mc%d.mat',  true;
};

methods_names = methods_dict(:,1);
methods_files = methods_dict(:,2);
methods_hastr = cell2mat(methods_dict(:,3));

ip_prefixes = {'IP-DAC','IP-AC'};
tp_prefixes = {'TP-DAC','TP-AC'};

% M-ablation file patterns (for Low/Medium/High rows)
agg_to_file = containers.Map( ...
    {'MOE','GPOE','POE','BCM','RBCM'}, ...
    {'moe','gpoe','poe','bcm','rbcm'});

%% ===================== Metrics to show =====================
% Each row: {display_name, field_in_mat, format_str, is_ip_only, is_tp_only}
metrics = {
    'SMSE',                  'smse',                '%.4f',  false, false;
    'RMSE',                  'rmse',                '%.4f',  false, false;
    'MSLL',                  'msll',                '%.4f',  false, false;
    'Train Time (ms/pt)',    't_train_per_point',   '%.4f',  false, false;
    'Test Time (ms/pt)',     't_test_per_point',    '%.4f',  false, false;
    'Trigger Times',         'comm_train',          '%.1f',  true,  false;
    'Consensus Iter',        'comm_test',           '%.0f',  false, true;
    'Trigger Ratio (%)',     'trigger_ratio_train', '%.1f',  true,  false;
};
n_metrics = size(metrics, 1);

%% ===================== Load main-experiment results =====================
% For non-IP methods, read from main experiment files (no M suffix)
fprintf('Loading main experiment results...\n');

n_methods = numel(methods_names);
n_ds      = numel(datasets);

mean_main = NaN(n_ds, n_methods, n_metrics);
std_main  = NaN(n_ds, n_methods, n_metrics);

for mi = 1:n_methods
    is_ip_method_main = startsWith(methods_names{mi}, 'IP-DAC') || ...
                        startsWith(methods_names{mi}, 'IP-AC');
    is_tp_method_main = startsWith(methods_names{mi}, 'TP-DAC') || ...
                        startsWith(methods_names{mi}, 'TP-AC');

    % determine prefix for this method
    parts = strsplit(methods_names{mi}, '-');
    if numel(parts) >= 3
        prefix = [parts{1} '-' parts{2}];
    else
        prefix = parts{1};
    end

    for di = 1:n_ds
        dname = datasets{di};
        folder = fullfile(ProjectRoot, 'Result', 'Dataset', dname);
        vals = NaN(n_metrics, n_mc);

        for ci = 1:n_mc
            mc = mc_list(ci);
            if methods_hastr(mi)
                fname = sprintf(methods_files{mi}, tr_tag, mc);
            else
                fname = sprintf(methods_files{mi}, mc);
            end
            fpath = fullfile(folder, fname);
            if ~isfile(fpath), continue; end
            try
                S = load(fpath);
                for ki = 1:n_metrics
                    field = metrics{ki, 2};
                    is_ip_only = metrics{ki, 4};
                    is_tp_only = metrics{ki, 5};
                    if is_ip_only && ~is_ip_method_main, continue; end
                    if is_tp_only && ~is_tp_method_main, continue; end
                    if isfield(S, field)
                        v = S.(field);
                        if isnumeric(v) && isscalar(v)
                            % Trigger ratio stored as fraction, convert to %
                            if strcmp(field, 'trigger_ratio_train')
                                v = v * 100;
                            end
                            vals(ki, ci) = v;
                        end
                    end
                end
            catch
            end
        end

        for ki = 1:n_metrics
            mean_main(di, mi, ki) = mean(vals(ki,:), 'omitnan');
            std_main(di, mi, ki)  = std(vals(ki,:),  'omitnan');
        end
    end
end

%% ===================== Generate HTML =====================

OutPath = fullfile(ProjectRoot, 'Result', 'Dataset', 'Summary_Table.html');
fid = fopen(OutPath, 'w', 'n', 'UTF-8');

fprintf(fid, '<!DOCTYPE html>\n<html>\n<head>\n');
fprintf(fid, '<meta charset="UTF-8">\n');
fprintf(fid, '<title>GP Methods Summary</title>\n');
fprintf(fid, '<style>\n');
fprintf(fid, 'body{font-family:Arial,sans-serif;font-size:13px;margin:20px;}\n');
fprintf(fid, 'h2{color:#1F4E79;font-size:15px;margin-top:30px;}\n');
fprintf(fid, 'table{border-collapse:collapse;width:100%%;margin-bottom:30px;}\n');
fprintf(fid, 'th,td{border:1px solid #BDD7EE;padding:3px 5px;text-align:center;white-space:nowrap;}\n');
fprintf(fid, 'th{background:#BDD7EE;font-weight:bold;font-size:13px;}\n');
fprintf(fid, 'th.agg{background:#1F4E79;color:white;}\n');
fprintf(fid, 'th.method{background:#2E75B6;color:white;}\n');
fprintf(fid, 'td.agg{background:#F2F2F2;font-weight:bold;text-align:center;}\n');
fprintf(fid, 'td.method{font-weight:bold;background:#FAFAFA;text-align:center;}\n');
fprintf(fid, 'td.mlevel{font-size:11px;color:#555;text-align:center;}\n');
fprintf(fid, 'td.val{font-size:13px;text-align:center;}\n');
fprintf(fid, 'td.dash{color:#AAA;text-align:center;}\n');
fprintf(fid, 'td.best{font-size:13px;text-align:center;font-weight:bold;background:#FFF2CC;color:#7F6000;}\n');
fprintf(fid, '</style>\n</head>\n<body>\n');

fprintf(fid, '<h1 style="color:#1F4E79;">GP Methods Summary &mdash; Train=40%%, MC=%d</h1>\n', n_mc);

for ki = 1:n_metrics
    metric_name  = metrics{ki, 1};
    metric_field = metrics{ki, 2};
    fmt          = metrics{ki, 3};
    is_ip_only   = metrics{ki, 4};
    is_tp_only   = metrics{ki, 5};
    is_ratio     = strcmp(metric_field, 'trigger_ratio_train');

    fprintf(fid, '<h2>%s &nbsp; (Train=40%%, MC=%d)</h2>\n', metric_name, n_mc);
    fprintf(fid, '<table>\n');

    % ---- Header row 1: Aggregation | Method | M Level | datasets ----
    fprintf(fid, '<tr>');
    fprintf(fid, '<th class="agg" rowspan="1">Aggregation</th>');
    fprintf(fid, '<th class="method" rowspan="1">Method</th>');
    fprintf(fid, '<th class="method" rowspan="1">M Level</th>');
    for di = 1:n_ds
        fprintf(fid, '<th colspan="1">%s</th>', datasets{di});
    end
    fprintf(fid, '</tr>\n');

    for a_idx = 1:numel(agg_list)
        agg       = agg_list{a_idx};
        agg_lower = agg_to_file(agg);

        % Count total rows for this agg block (for rowspan)
        total_rows = 0;
        for col = 1:numel(col_prefixes)
            prefix = col_prefixes{col};
            is_ip  = ismember(prefix, ip_prefixes);
            if is_ip && ~is_ip_only && ~is_tp_only
                total_rows = total_rows + 3;  % Low/Medium/High
            else
                total_rows = total_rows + 1;
            end
        end

        % Pre-compute best M level per prefix (IP-DAC / IP-AC) per dataset
        % best_lv_per_prefix{prefix_idx}(di) = best lv index (1/2/3) for that prefix+dataset
        is_larger_better = strcmp(metric_field,'msll') || strcmp(metric_field,'trigger_ratio_train');
        best_lv_per_prefix = cell(numel(col_prefixes), 1);
        for col = 1:numel(col_prefixes)
            prefix_tmp = col_prefixes{col};
            if ~ismember(prefix_tmp, ip_prefixes)
                best_lv_per_prefix{col} = zeros(1, n_ds);  % 0 = no highlight
                continue;
            end
            % Get agg file key for this prefix
            best_lvs = zeros(1, n_ds);
            for di = 1:n_ds
                dname_tmp = datasets{di};
                best_v = NaN;
                best_lv = 0;
                for lv = 1:3
                    if strcmp(dname_tmp,'PUMADYN32NM')
                        M_tmp = LMH_puma(lv);
                    else
                        M_tmp = LMH_default(lv);
                    end
                    [mv_tmp,~] = load_mablation(ProjectRoot, dname_tmp, agg_lower, ...
                        M_tmp, tr_tag, n_mc, metric_field, is_ratio);
                    if ~isnan(mv_tmp)
                        if isnan(best_v) || ...
                           (is_larger_better && mv_tmp > best_v) || ...
                           (~is_larger_better && mv_tmp < best_v)
                            best_v  = mv_tmp;
                            best_lv = lv;
                        end
                    end
                end
                best_lvs(di) = best_lv;
            end
            best_lv_per_prefix{col} = best_lvs;
        end

        first_row = true;

        for col = 1:numel(col_prefixes)
            prefix      = col_prefixes{col};
            target_name = sprintf('%s-%s', prefix, agg);
            mi          = find(strcmp(methods_names, target_name), 1);
            is_ip       = ismember(prefix, ip_prefixes);
            is_tp       = ismember(prefix, tp_prefixes);

            % For IP methods with normal metrics: show Low/Medium/High rows
            show_lmh = is_ip && ~is_ip_only && ~is_tp_only && ~is_ratio;

            if show_lmh
                for lv = 1:3
                    fprintf(fid, '<tr>');

                    % Aggregation cell (only on very first row)
                    if first_row && lv == 1
                        fprintf(fid, '<td class="agg" rowspan="%d">%s</td>', total_rows, agg);
                        first_row = false;
                    end

                    % Method cell (only on first level row)
                    if lv == 1
                        fprintf(fid, '<td class="method" rowspan="3">%s</td>', prefix);
                    end

                    % M Level label
                    lmh_labels = {'Low','Medium','High'};
                    fprintf(fid, '<td class="mlevel">%s</td>', lmh_labels{lv});

                    % Data cells from M-ablation files
                    for di = 1:n_ds
                        dname = datasets{di};
                        if strcmp(dname, 'PUMADYN32NM')
                            M_val = LMH_puma(lv);
                        else
                            M_val = LMH_default(lv);
                        end

                        [mv, sv] = load_mablation(ProjectRoot, dname, agg_lower, ...
                            M_val, tr_tag, n_mc, metric_field, is_ratio);

                        if isnan(mv)
                            fprintf(fid, '<td class="dash">-</td>');
                        else
                            val_str = sprintf(fmt, mv);
                            std_str = sprintf(fmt, sv);
                            % Highlight only the single best M level for this prefix+dataset
                            is_best = (best_lv_per_prefix{col}(di) == lv);
                            if is_best
                                fprintf(fid, '<td class="best">%s&plusmn;%s</td>', val_str, std_str);
                            else
                                fprintf(fid, '<td class="val">%s&plusmn;%s</td>', val_str, std_str);
                            end
                        end
                    end
                    fprintf(fid, '</tr>\n');
                end

            else
                % Single row
                fprintf(fid, '<tr>');

                if first_row
                    fprintf(fid, '<td class="agg" rowspan="%d">%s</td>', total_rows, agg);
                    first_row = false;
                end

                fprintf(fid, '<td class="method">%s</td>', prefix);
                fprintf(fid, '<td class="mlevel">-</td>');

                for di = 1:n_ds
                    if isempty(mi)
                        fprintf(fid, '<td class="dash">-</td>');
                        continue;
                    end

                    % Hide trigger metrics for non-IP / non-TP
                    if is_ip_only && ~is_ip
                        fprintf(fid, '<td class="dash">-</td>');
                        continue;
                    end
                    if is_tp_only && ~is_tp
                        fprintf(fid, '<td class="dash">-</td>');
                        continue;
                    end

                    mv = mean_main(di, mi, ki);
                    sv = std_main(di, mi, ki);

                    if isnan(mv)
                        fprintf(fid, '<td class="dash">-</td>');
                    else
                        val_str = sprintf(fmt, mv);
                        std_str = sprintf(fmt, sv);
                        fprintf(fid, '<td class="val">%s&plusmn;%s</td>', val_str, std_str);
                    end
                end
                fprintf(fid, '</tr>\n');
            end
        end
    end

    fprintf(fid, '</table>\n\n');
end

fprintf(fid, '</body>\n</html>\n');
fclose(fid);

fprintf('\nDone! HTML saved to:\n%s\n', OutPath);
fprintf('Open in browser (Chrome/Edge), then copy-paste into Word.\n');

%% ===================== Local function (must be at end of script) =====================
function [mv, sv] = load_mablation(ProjectRoot, dname, agg_lower, M_val, tr_tag, n_mc, field, is_ratio)
    folder = fullfile(ProjectRoot, 'Result', 'Dataset', dname);
    vals = NaN(1, n_mc);
    for mc = 1:n_mc
        fname = sprintf('%s_M%d_tr%d_mc%d.mat', agg_lower, M_val, tr_tag, mc);
        fpath = fullfile(folder, fname);
        if ~isfile(fpath), continue; end
        try
            S = load(fpath);
            if isfield(S, field)
                v = S.(field);
                if isnumeric(v) && isscalar(v)
                    if is_ratio, v = v * 100; end
                    vals(mc) = v;
                end
            end
        catch
        end
    end
    mv = mean(vals, 'omitnan');
    sv = std(vals,  'omitnan');
end