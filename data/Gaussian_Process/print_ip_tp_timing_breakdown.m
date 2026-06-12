function print_ip_tp_timing_breakdown()
% PRINT_IP_TP_TIMING_BREAKDOWN_MSPT
% Reads saved .mat files and prints two timing-breakdown tables:
%   1) IP-DAC/IP-AC training-side breakdown, unit = ms / training point
%   2) TP-DAC/TP-AC test-side breakdown, unit = ms / test point

clc;

%datasets = {'KIN40K', 'POL', 'PUMADYN32NM', 'SARCOS'};
datasets = { 'SARCOS'};
aggs     = {'moe', 'gpoe', 'poe', 'bcm', 'rbcm'};

train_ratio = 0.4;
tr_tag = round(train_ratio * 100);
mc_seeds = 1:10;

ResultRoot = fullfile('Result', 'Dataset');

%% ========================================================================
% Table 1: IP training-side timing, normalized by N_train
% ========================================================================
fprintf('\n============================================================================================================================\n');
fprintf('Training-side timing breakdown for IP-DAC / IP-AC  (Train=%d%%, MC=%d)\n', tr_tag, numel(mc_seeds));
fprintf('Unit: ms per training point\n');
fprintf('============================================================================================================================\n');
fprintf('%-6s %-13s %-8s %16s %20s %14s %18s %14s\n', ...
    'Agg', 'Dataset', 'Mode', 'LocalGP train', 'Inducing local pred.', ...
    'Consensus', 'MaskedGP train', 'Total train');
fprintf('----------------------------------------------------------------------------------------------------------------------------\n');

for ai = 1:numel(aggs)
    agg = aggs{ai};

    for di = 1:numel(datasets)
        DatasetName = datasets{di};
        ResultFolder = fullfile(ResultRoot, DatasetName);

        ip_dac_stats = collect_ip_stats_mspt(ResultFolder, agg, '', tr_tag, mc_seeds);
        print_ip_row(agg, DatasetName, 'IP-DAC', ip_dac_stats);

        ip_ac_stats = collect_ip_stats_mspt(ResultFolder, agg, '_ac', tr_tag, mc_seeds);
        print_ip_row(agg, DatasetName, 'IP-AC', ip_ac_stats);

        fprintf('----------------------------------------------------------------------------------------------------------------------------\n');
    end

    fprintf('============================================================================================================================\n');
end


%% ========================================================================
% Table 2: TP test-side timing, normalized by N_eval
% ========================================================================
fprintf('\n==============================================================================================================\n');
fprintf('Test-side timing breakdown for TP-DAC / TP-AC  (Train=%d%%, MC=%d)\n', tr_tag, numel(mc_seeds));
fprintf('Unit: ms per test point\n');
fprintf('==============================================================================================================\n');
fprintf('%-6s %-13s %-8s %22s %14s %14s\n', ...
    'Agg', 'Dataset', 'Mode', 'Test local prediction', 'Consensus', 'Total test');
fprintf('--------------------------------------------------------------------------------------------------------------\n');

for ai = 1:numel(aggs)
    agg = aggs{ai};

    for di = 1:numel(datasets)
        DatasetName = datasets{di};
        ResultFolder = fullfile(ResultRoot, DatasetName);

        tp_dac_stats = collect_tp_stats_mspt(ResultFolder, agg, '_tp', tr_tag, mc_seeds);
        print_tp_row(agg, DatasetName, 'TP-DAC', tp_dac_stats);

        tp_ac_stats = collect_tp_stats_mspt_multi_suffix(ResultFolder, agg, {'_ac_tp','_tp_ac'}, tr_tag, mc_seeds);
        print_tp_row(agg, DatasetName, 'TP-AC', tp_ac_stats);

        fprintf('--------------------------------------------------------------------------------------------------------------\n');
    end

    fprintf('==============================================================================================================\n');
end

fprintf('\nNote:\n');
fprintf('  IP table: all entries are normalized by N_train and represent training-side cost.\n');
fprintf('  TP table: entries are normalized by N_eval and represent test-side aggregation cost.\n');
fprintf('  LocalGP training is shared preprocessing and is not TP-specific test-time cost.\n\n');

end


function stats = collect_ip_stats_mspt(ResultFolder, agg, suffix, tr_tag, mc_seeds)

local_gp   = [];
local_pred = [];
consensus  = [];
maskedgp   = [];
total      = [];

for seed = mc_seeds
    file = fullfile(ResultFolder, sprintf('%s%s_tr%d_mc%d.mat', agg, suffix, tr_tag, seed));
    if ~exist(file, 'file'), continue; end

    d = load(file);
    N_train = infer_N_train(d);

    if isnan(N_train) || N_train <= 0
        warning('Cannot infer N_train for file: %s. Skipping.', file);
        continue;
    end

    scale = 1000 / N_train;   % seconds -> ms/train point

    local_gp(end+1)   = getfield_default(d, 't_ip_local_gp_train', NaN)      * scale; %#ok<AGROW>
    local_pred(end+1) = getfield_default(d, 't_ip_inducing_prediction', NaN) * scale; %#ok<AGROW>
    consensus(end+1)  = getfield_default(d, 't_ip_consensus', ...
                         getfield_default(d, 't_consensus', NaN))            * scale; %#ok<AGROW>
    maskedgp(end+1)   = getfield_default(d, 't_ip_maskedgp_train', NaN)      * scale; %#ok<AGROW>
    total(end+1)      = getfield_default(d, 't_ip_total_train', ...
                         getfield_default(d, 't_train_total', ...
                         getfield_default(d, 't_total_train', NaN)))         * scale; %#ok<AGROW>
end

stats.local_gp   = summarize(local_gp);
stats.local_pred = summarize(local_pred);
stats.consensus  = summarize(consensus);
stats.maskedgp   = summarize(maskedgp);
stats.total      = summarize(total);

end


function stats = collect_tp_stats_mspt(ResultFolder, agg, suffix, tr_tag, mc_seeds)

local_pred = [];
consensus  = [];
total      = [];

for seed = mc_seeds
    file = fullfile(ResultFolder, sprintf('%s%s_tr%d_mc%d.mat', agg, suffix, tr_tag, seed));
    if ~exist(file, 'file'), continue; end

    d = load(file);
    N_eval = infer_N_eval(d);

    if isnan(N_eval) || N_eval <= 0
        warning('Cannot infer N_eval for file: %s. Skipping.', file);
        continue;
    end

    scale = 1000 / N_eval;   % seconds -> ms/test point

    local_pred(end+1) = getfield_default(d, 't_tp_test_local_prediction', NaN) * scale; %#ok<AGROW>
    consensus(end+1)  = getfield_default(d, 't_tp_consensus', ...
                         getfield_default(d, 't_consensus', NaN))              * scale; %#ok<AGROW>
    total(end+1)      = getfield_default(d, 't_tp_total_test', ...
                         getfield_default(d, 't_test_total', NaN))             * scale; %#ok<AGROW>
end

stats.local_pred = summarize(local_pred);
stats.consensus  = summarize(consensus);
stats.total      = summarize(total);

end


function stats = collect_tp_stats_mspt_multi_suffix(ResultFolder, agg, suffix_list, tr_tag, mc_seeds)

local_pred = [];
consensus  = [];
total      = [];

for seed = mc_seeds
    found = false;

    for si = 1:numel(suffix_list)
        suffix = suffix_list{si};
        file = fullfile(ResultFolder, sprintf('%s%s_tr%d_mc%d.mat', agg, suffix, tr_tag, seed));
        if exist(file, 'file')
            found = true;
            break;
        end
    end

    if ~found, continue; end

    d = load(file);
    N_eval = infer_N_eval(d);

    if isnan(N_eval) || N_eval <= 0
        warning('Cannot infer N_eval for file: %s. Skipping.', file);
        continue;
    end

    scale = 1000 / N_eval;

    local_pred(end+1) = getfield_default(d, 't_tp_test_local_prediction', NaN) * scale; %#ok<AGROW>
    consensus(end+1)  = getfield_default(d, 't_tp_consensus', ...
                         getfield_default(d, 't_consensus', NaN))              * scale; %#ok<AGROW>
    total(end+1)      = getfield_default(d, 't_tp_total_test', ...
                         getfield_default(d, 't_test_total', NaN))             * scale; %#ok<AGROW>
end

stats.local_pred = summarize(local_pred);
stats.consensus  = summarize(consensus);
stats.total      = summarize(total);

end


function print_ip_row(agg, DatasetName, mode, stats)

fprintf('%-6s %-13s %-8s %16s %20s %14s %18s %14s\n', ...
    upper(agg), DatasetName, mode, ...
    format_summary(stats.local_gp), ...
    format_summary(stats.local_pred), ...
    format_summary(stats.consensus), ...
    format_summary(stats.maskedgp), ...
    format_summary(stats.total));

end


function print_tp_row(agg, DatasetName, mode, stats)

fprintf('%-6s %-13s %-8s %22s %14s %14s\n', ...
    upper(agg), DatasetName, mode, ...
    format_summary(stats.local_pred), ...
    format_summary(stats.consensus), ...
    format_summary(stats.total));

end


function N_train = infer_N_train(d)

N_train = NaN;

if isfield(d, 'N_train')
    N_train = d.N_train;
    return;
end

if isfield(d, 't_train_total') && isfield(d, 't_train_per_point') && d.t_train_per_point > 0
    N_train = d.t_train_total * 1000 / d.t_train_per_point;
    return;
end

if isfield(d, 't_ip_total_train') && isfield(d, 't_train_per_point') && d.t_train_per_point > 0
    N_train = d.t_ip_total_train * 1000 / d.t_train_per_point;
    return;
end

end


function N_eval = infer_N_eval(d)

N_eval = NaN;

if isfield(d, 'N_eval')
    N_eval = d.N_eval;
    return;
end

if isfield(d, 't_test_total') && isfield(d, 't_test_per_point') && d.t_test_per_point > 0
    N_eval = d.t_test_total * 1000 / d.t_test_per_point;
    return;
end

if isfield(d, 't_tp_total_test') && isfield(d, 't_test_per_point') && d.t_test_per_point > 0
    N_eval = d.t_tp_total_test * 1000 / d.t_test_per_point;
    return;
end

end


function s = summarize(x)

x = x(~isnan(x));

if isempty(x)
    s.mean = NaN;
    s.std  = NaN;
    s.n    = 0;
else
    s.mean = mean(x);
    if numel(x) <= 1
        s.std = 0;
    else
        s.std = std(x);
    end
    s.n = numel(x);
end

end


function str = format_summary(s)

if s.n == 0 || isnan(s.mean)
    str = '-';
else
    str = sprintf('%.3f±%.3f', s.mean, s.std);
end

end


function value = getfield_default(s, field_name, default_value)

if isfield(s, field_name)
    value = s.(field_name);
else
    value = default_value;
end

end
