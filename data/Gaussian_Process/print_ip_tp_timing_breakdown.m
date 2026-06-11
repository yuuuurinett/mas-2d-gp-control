function print_ip_tp_timing_breakdown()
% PRINT_IP_TP_TIMING_BREAKDOWN
% Print timing breakdown table for IP-DAC / IP-AC / TP-DAC / TP-AC.
%
% This script only reads saved .mat result files. It does not rerun experiments.
%
% Table columns:
%   Local prediction:
%       IP methods: t_ip_inducing_prediction
%       TP methods: t_tp_test_local_prediction
%
%   Consensus:
%       IP methods: t_ip_consensus
%       TP methods: t_tp_consensus
%
%   MaskedGP train:
%       IP methods: t_ip_maskedgp_train
%       TP methods: shown as '-'
%
%   Total:
%       IP methods: t_ip_total_train
%       TP methods: t_tp_total_test
%
% Units: seconds, averaged over MC seeds.

clc;

%% ===================== User settings =====================
datasets = {'KIN40K', 'POL', 'PUMADYN32NM', 'SARCOS'};
aggs     = {'moe', 'gpoe', 'poe', 'bcm', 'rbcm'};

train_ratio = 0.4;
tr_tag = round(train_ratio * 100);

% Change this if you use a different MC range.
mc_seeds = 1:10;

% Assumes this script is run from Gaussian_Process folder.
ResultRoot = fullfile('Result', 'Dataset');

%% ===================== Print table =====================
fprintf('\n============================================================================================================================\n');
fprintf('Timing Breakdown for IP-DAC / IP-AC / TP-DAC / TP-AC  (Train=%d%%, MC=%d)\n', ...
    tr_tag, numel(mc_seeds));
fprintf('============================================================================================================================\n');

fprintf('%-6s %-13s %-8s %18s %14s %18s %14s\n', ...
    'Agg', 'Dataset', 'Mode', 'Local prediction', 'Consensus', 'MaskedGP train', 'Total');
fprintf('----------------------------------------------------------------------------------------------------------------------------\n');

for ai = 1:numel(aggs)
    agg = aggs{ai};

    for di = 1:numel(datasets)
        DatasetName = datasets{di};
        ResultFolder = fullfile(ResultRoot, DatasetName);

        % IP-DAC file name: poe_tr40_mc1.mat
        ip_dac_stats = collect_ip_stats(ResultFolder, agg, '', tr_tag, mc_seeds);
        print_row(agg, DatasetName, 'IP-DAC', ip_dac_stats);

        % IP-AC file name: poe_ac_tr40_mc1.mat
        ip_ac_stats = collect_ip_stats(ResultFolder, agg, '_ac', tr_tag, mc_seeds);
        print_row(agg, DatasetName, 'IP-AC', ip_ac_stats);

        % TP-DAC file name: poe_tp_tr40_mc1.mat
        tp_dac_stats = collect_tp_stats(ResultFolder, agg, '_tp', tr_tag, mc_seeds);
        print_row(agg, DatasetName, 'TP-DAC', tp_dac_stats);

        % TP-AC file name: poe_tp_ac_tr40_mc1.mat
        tp_ac_stats = collect_tp_stats(ResultFolder, agg, '_tp_ac', tr_tag, mc_seeds);
        print_row(agg, DatasetName, 'TP-AC', tp_ac_stats);

        fprintf('----------------------------------------------------------------------------------------------------------------------------\n');
    end

    fprintf('============================================================================================================================\n');
end

fprintf('\nNote:\n');
fprintf('  IP-DAC/IP-AC Total = training-side total time.\n');
fprintf('  TP-DAC/TP-AC Total = test-side total time.\n');
fprintf('  MaskedGP train applies only to IP methods; TP methods show ''.''.\n');
fprintf('  All values are in seconds and reported as mean±std over MC seeds.\n\n');

end


%% ========================================================================
% Collect IP-DAC / IP-AC timing
% ========================================================================
function stats = collect_ip_stats(ResultFolder, agg, suffix, tr_tag, mc_seeds)
% suffix:
%   ''    for IP-DAC file: poe_tr40_mc1.mat
%   '_ac' for IP-AC  file: poe_ac_tr40_mc1.mat

local_pred = [];
consensus  = [];
maskedgp   = [];
total      = [];

for seed = mc_seeds
    file = fullfile(ResultFolder, sprintf('%s%s_tr%d_mc%d.mat', agg, suffix, tr_tag, seed));

    if ~exist(file, 'file')
        continue;
    end

    d = load(file);

    % IP local prediction at inducing points.
    t_local_pred = getfield_default(d, 't_ip_inducing_prediction', NaN);

    % IP consensus stage.
    t_consensus = getfield_default(d, 't_ip_consensus', ...
                  getfield_default(d, 't_consensus', NaN));

    % IP MaskedGP training only.
    t_maskedgp = getfield_default(d, 't_ip_maskedgp_train', NaN);

    % IP training-side total only.
    t_total = getfield_default(d, 't_ip_total_train', ...
              getfield_default(d, 't_total_train', ...
              getfield_default(d, 't_train_total', NaN)));

    local_pred(end+1) = t_local_pred; %#ok<AGROW>
    consensus(end+1)  = t_consensus;  %#ok<AGROW>
    maskedgp(end+1)   = t_maskedgp;   %#ok<AGROW>
    total(end+1)      = t_total;      %#ok<AGROW>
end

stats.local_pred = summarize(local_pred);
stats.consensus  = summarize(consensus);
stats.maskedgp   = summarize(maskedgp);
stats.total      = summarize(total);
stats.has_maskedgp = true;

end


%% ========================================================================
% Collect TP-DAC / TP-AC timing
% ========================================================================
function stats = collect_tp_stats(ResultFolder, agg, suffix, tr_tag, mc_seeds)
% suffix:
%   '_tp'    for TP-DAC file: poe_tp_tr40_mc1.mat
%   '_tp_ac' for TP-AC  file: poe_tp_ac_tr40_mc1.mat

local_pred = [];
consensus  = [];
total      = [];

for seed = mc_seeds
    file = fullfile(ResultFolder, sprintf('%s%s_tr%d_mc%d.mat', agg, suffix, tr_tag, seed));

    if ~exist(file, 'file')
        continue;
    end

    d = load(file);

    % TP local prediction at test points.
    t_local_pred = getfield_default(d, 't_tp_test_local_prediction', NaN);

    % TP consensus stage.
    t_consensus = getfield_default(d, 't_tp_consensus', ...
                  getfield_default(d, 't_consensus', NaN));

    % TP test-side total.
    t_total = getfield_default(d, 't_tp_total_test', ...
              getfield_default(d, 't_test_total', NaN));

    local_pred(end+1) = t_local_pred; %#ok<AGROW>
    consensus(end+1)  = t_consensus;  %#ok<AGROW>
    total(end+1)      = t_total;      %#ok<AGROW>
end

stats.local_pred = summarize(local_pred);
stats.consensus  = summarize(consensus);
stats.maskedgp   = [];       % TP methods do not train MaskedGP.
stats.total      = summarize(total);
stats.has_maskedgp = false;

end


%% ========================================================================
% Print one row
% ========================================================================
function print_row(agg, DatasetName, mode, stats)

local_str = format_summary(stats.local_pred);
cons_str  = format_summary(stats.consensus);

if stats.has_maskedgp
    maskedgp_str = format_summary(stats.maskedgp);
else
    maskedgp_str = '-';
end

total_str = format_summary(stats.total);

fprintf('%-6s %-13s %-8s %18s %14s %18s %14s\n', ...
    upper(agg), DatasetName, mode, local_str, cons_str, maskedgp_str, total_str);

end


%% ========================================================================
% Summary helper
% ========================================================================
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
% Format helper
% ========================================================================
function str = format_summary(s)

if s.n == 0 || isnan(s.mean)
    str = '-';
else
    str = sprintf('%.3f%c%.3f', s.mean, char(177), s.std);
end

end


%% ========================================================================
% Safe field getter
% ========================================================================
function value = getfield_default(s, field_name, default_value)

if isfield(s, field_name)
    value = s.(field_name);
else
    value = default_value;
end

end
