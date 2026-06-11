function print_ip_tp_timing_breakdown()
% PRINT_IP_TP_TIMING_BREAKDOWN
% Print timing breakdown table for IP-DAC / IP-AC / TP-DAC / TP-AC.
%
% The script only reads saved .mat result files. It does not rerun experiments.
%
% Table columns:
%   LocalGP train      IP: t_ip_local_gp_train, TP: t_tp_local_gp_train
%   Local prediction   IP: t_ip_inducing_prediction, TP: t_tp_test_local_prediction
%   Consensus          IP: t_ip_consensus, TP: t_tp_consensus
%   MaskedGP train     IP: t_ip_maskedgp_train, TP: '-'
%   Total              IP: t_ip_total_train, TP: t_tp_total_test
%
% Put this file under Gaussian_Process and run:
%   print_ip_tp_timing_breakdown

clc;

%% Settings
%datasets = {'KIN40K', 'POL', 'PUMADYN32NM', 'SARCOS'};
datasets = { 'SARCOS'};
aggs     = {'moe', 'gpoe', 'poe', 'bcm', 'rbcm'};

train_ratio = 0.4;
tr_tag = round(train_ratio * 100);
mc_seeds = 1;
%mc_seeds = 1:10;


ResultRoot = fullfile('Result', 'Dataset');

fprintf('\n====================================================================================================================================\n');
fprintf('Timing Breakdown for IP-DAC / IP-AC / TP-DAC / TP-AC  (Train=%d%%, MC=%d)\n', tr_tag, numel(mc_seeds));
fprintf('====================================================================================================================================\n');
fprintf('%-6s %-13s %-8s %16s %18s %14s %18s %14s\n', ...
    'Agg', 'Dataset', 'Mode', 'LocalGP train', 'Local prediction', ...
    'Consensus', 'MaskedGP train', 'Total');
fprintf('------------------------------------------------------------------------------------------------------------------------------------\n');

for ai = 1:numel(aggs)
    agg = aggs{ai};

    for di = 1:numel(datasets)
        DatasetName = datasets{di};
        ResultFolder = fullfile(ResultRoot, DatasetName);

        % IP-DAC: e.g. poe_tr40_mc1.mat
        ip_dac_stats = collect_ip_stats(ResultFolder, agg, '', tr_tag, mc_seeds);
        print_row(agg, DatasetName, 'IP-DAC', ip_dac_stats);

        % IP-AC: e.g. poe_ac_tr40_mc1.mat
        ip_ac_stats = collect_ip_stats(ResultFolder, agg, '_ac', tr_tag, mc_seeds);
        print_row(agg, DatasetName, 'IP-AC', ip_ac_stats);

        % TP-DAC: try possible naming conventions
        tp_dac_stats = collect_tp_stats_multi_suffix(ResultFolder, agg, {'_tp', '_testpoint'}, tr_tag, mc_seeds);
        print_row(agg, DatasetName, 'TP-DAC', tp_dac_stats);

        % TP-AC: try possible naming conventions
        tp_ac_stats = collect_tp_stats_multi_suffix(ResultFolder, agg, {'_tp_ac', '_ac_tp', '_testpoint_ac', '_ac_testpoint'}, tr_tag, mc_seeds);
        print_row(agg, DatasetName, 'TP-AC', tp_ac_stats);

        fprintf('------------------------------------------------------------------------------------------------------------------------------------\n');
    end

    fprintf('====================================================================================================================================\n');
end

fprintf('\nNotes:\n');
fprintf('  IP-DAC/IP-AC Total = training-side total time.\n');
fprintf('  TP-DAC/TP-AC Total = test-side total time.\n');
fprintf('  MaskedGP train applies only to IP methods.\n');
fprintf('  TP-AC file names are detected using multiple possible suffixes.\n\n');

end

%% ========================================================================
function stats = collect_ip_stats(ResultFolder, agg, suffix, tr_tag, mc_seeds)
% suffix:
%   ''    for IP-DAC file: poe_tr40_mc1.mat
%   '_ac' for IP-AC  file: poe_ac_tr40_mc1.mat

local_gp   = [];
local_pred = [];
consensus  = [];
gp_step    = [];
total      = [];

for seed = mc_seeds
    file = fullfile(ResultFolder, sprintf('%s%s_tr%d_mc%d.mat', agg, suffix, tr_tag, seed));

    if ~exist(file, 'file')
        continue;
    end

    d = load(file);

    t_local_gp = getfield_default(d, 't_ip_local_gp_train', NaN);
    t_local_pred = getfield_default(d, 't_ip_inducing_prediction', NaN);
    t_consensus = getfield_default(d, 't_ip_consensus', getfield_default(d, 't_consensus', NaN));
    t_gp_step = getfield_default(d, 't_ip_maskedgp_train', NaN);
    t_total = getfield_default(d, 't_ip_total_train', ...
              getfield_default(d, 't_total_train', ...
              getfield_default(d, 't_train_total', NaN)));

    local_gp(end+1)   = t_local_gp;    %#ok<AGROW>
    local_pred(end+1) = t_local_pred;  %#ok<AGROW>
    consensus(end+1)  = t_consensus;   %#ok<AGROW>
    gp_step(end+1)    = t_gp_step;     %#ok<AGROW>
    total(end+1)      = t_total;       %#ok<AGROW>
end

stats.local_gp   = summarize(local_gp);
stats.local_pred = summarize(local_pred);
stats.consensus  = summarize(consensus);
stats.gp_step    = summarize(gp_step);
stats.total      = summarize(total);
stats.has_gp_step = true;

end

%% ========================================================================
function stats = collect_tp_stats_multi_suffix(ResultFolder, agg, suffix_list, tr_tag, mc_seeds)
% Try several file suffix conventions and collect whichever exists.

local_gp   = [];
local_pred = [];
consensus  = [];
total      = [];

for seed = mc_seeds
    file = '';
    for si = 1:numel(suffix_list)
        candidate = fullfile(ResultFolder, sprintf('%s%s_tr%d_mc%d.mat', agg, suffix_list{si}, tr_tag, seed));
        if exist(candidate, 'file')
            file = candidate;
            break;
        end
    end

    if isempty(file)
        continue;
    end

    d = load(file);

    t_local_gp = getfield_default(d, 't_tp_local_gp_train', NaN);
    t_local_pred = getfield_default(d, 't_tp_test_local_prediction', NaN);
    t_consensus = getfield_default(d, 't_tp_consensus', getfield_default(d, 't_consensus', NaN));
    t_total = getfield_default(d, 't_tp_total_test', getfield_default(d, 't_test_total', NaN));

    local_gp(end+1)   = t_local_gp;    %#ok<AGROW>
    local_pred(end+1) = t_local_pred;  %#ok<AGROW>
    consensus(end+1)  = t_consensus;   %#ok<AGROW>
    total(end+1)      = t_total;       %#ok<AGROW>
end

stats.local_gp   = summarize(local_gp);
stats.local_pred = summarize(local_pred);
stats.consensus  = summarize(consensus);
stats.gp_step    = summarize([]);
stats.total      = summarize(total);
stats.has_gp_step = false;

end

%% ========================================================================
function print_row(agg, DatasetName, mode, stats)

local_gp_str = format_summary(stats.local_gp);
local_str    = format_summary(stats.local_pred);
cons_str     = format_summary(stats.consensus);

if isfield(stats, 'has_gp_step') && stats.has_gp_step
    gp_str = format_summary(stats.gp_step);
else
    gp_str = '-';
end

total_str = format_summary(stats.total);

fprintf('%-6s %-13s %-8s %16s %18s %14s %18s %14s\n', ...
    upper(agg), DatasetName, mode, local_gp_str, local_str, ...
    cons_str, gp_str, total_str);

end

%% ========================================================================
function s = summarize(x)

x = x(~isnan(x));

if isempty(x)
    s.mean = NaN;
    s.std  = NaN;
    s.n    = 0;
else
    s.mean = mean(x);
    if numel(x) == 1
        s.std = 0;
    else
        s.std = std(x);
    end
    s.n = numel(x);
end

end

%% ========================================================================
function str = format_summary(s)

if s.n == 0 || isnan(s.mean)
    str = '-';
else
    str = sprintf('%.3f±%.3f', s.mean, s.std);
end

end

%% ========================================================================
function value = getfield_default(s, field_name, default_value)

if isfield(s, field_name)
    value = s.(field_name);
else
    value = default_value;
end

end
