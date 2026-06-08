%% plot_consensus.m
%   (1) IP-DAC tracking error to the average reference input
%   (2) IP-DAC inter-agent disagreement
%   (3) IP-AC error to the initial average value
%

clc; close all;

DatasetName = 'KIN40K';
base_method = 'poe';
train_ratio = 0.4;
seed = 1;

tr_tag = round(train_ratio * 100);
SaveFolder = fullfile('Result','Dataset',DatasetName);

dac_file = fullfile(SaveFolder, sprintf('%s_tr%d_mc%d.mat', base_method, tr_tag, seed));
ac_file  = fullfile(SaveFolder, sprintf('%s_ac_tr%d_mc%d.mat', base_method, tr_tag, seed));

if ~exist(dac_file,'file')
    error('DAC result file not found: %s\nRun run_inducingpoint_dataset(''%s'',''%s'',%.2f,%d) first.', ...
        dac_file, DatasetName, base_method, train_ratio, seed);
end
if ~exist(ac_file,'file')
    error('AC result file not found: %s\nRun run_inducingpoint_dataset(''%s'',''%s_ac'',%.2f,%d) first.', ...
        ac_file, DatasetName, base_method, train_ratio, seed);
end

d_dac = load(dac_file);
d_ac  = load(ac_file);

% Backward-compatible fallbacks in case an older result file is loaded.
if isfield(d_dac,'conv_dac_tracking_curve')
    dac_tracking = d_dac.conv_dac_tracking_curve;
else
    warning('conv_dac_tracking_curve not found. Using conv_curve_dac as fallback.');
    dac_tracking = d_dac.conv_curve_dac;
end

if isfield(d_dac,'conv_dac_disagreement_curve')
    dac_disagreement = d_dac.conv_dac_disagreement_curve;
elseif isfield(d_dac,'conv_curve_dac')
    dac_disagreement = d_dac.conv_curve_dac;
else
    error('No DAC convergence curve found in %s', dac_file);
end

if isfield(d_ac,'conv_ac_avg_error_curve')
    ac_avg_error = d_ac.conv_ac_avg_error_curve;
elseif isfield(d_ac,'conv_curve_ac')
    warning('conv_ac_avg_error_curve not found. Using conv_curve_ac as fallback.');
    ac_avg_error = d_ac.conv_curve_ac;
else
    error('No AC convergence curve found in %s', ac_file);
end

%% Figure 1: three scalar convergence-error tests
figure('Color','w','Position',[100,100,950,720]);

subplot(3,1,1);
semilogy(1:numel(dac_tracking), dac_tracking, 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Mean max error');
title('IP-DAC test 1: tracking error to average reference input');
grid on;

subplot(3,1,2);
semilogy(1:numel(dac_disagreement), dac_disagreement, 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Mean max disagreement');
title('IP-DAC test 2: inter-agent disagreement');
grid on;

subplot(3,1,3);
semilogy(1:numel(ac_avg_error), ac_avg_error, 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Mean max error');
title('IP-AC test: error to initial average');
grid on;

sgtitle(sprintf('Consensus diagnostics [%s, %s, train=%d%%, seed=%d]', ...
    DatasetName, upper(base_method), tr_tag, seed));

%% Figure 2: representative state trajectories for one dimension and one inducing point
% This makes the convergence visually intuitive: six agent curves should
% approach the reference line.

if isfield(d_dac,'dac_state_hist') && isfield(d_dac,'dac_ref_hist') && ~isempty(d_dac.dac_state_hist)
    figure('Color','w','Position',[120,120,900,360]);
    plot(d_dac.dac_state_hist, 'LineWidth', 1.2); hold on;
    plot(d_dac.dac_ref_hist, 'k--', 'LineWidth', 1.6);
    xlabel('Iteration');
    ylabel('Consensus state');
    title(sprintf('IP-DAC representative state: dimension r=%d, inducing point m=%d', ...
        d_dac.conv_r, d_dac.conv_m));
    legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','Average reference', ...
        'Location','bestoutside');
    grid on;
end

if isfield(d_ac,'ac_state_hist') && isfield(d_ac,'ac_ref_hist') && ~isempty(d_ac.ac_state_hist)
    figure('Color','w','Position',[140,140,900,360]);
    plot(d_ac.ac_state_hist, 'LineWidth', 1.2); hold on;
    plot(d_ac.ac_ref_hist, 'k--', 'LineWidth', 1.6);
    xlabel('Iteration');
    ylabel('Consensus state');
    title(sprintf('IP-AC representative state: dimension r=%d, inducing point m=%d', ...
        d_ac.conv_r, d_ac.conv_m));
    legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','Initial average', ...
        'Location','bestoutside');
    grid on;
end

%% Print a compact explanation
fprintf('\nLoaded files:\n  DAC: %s\n  AC : %s\n', dac_file, ac_file);
if isfield(d_dac,'comm_train')
    fprintf('DAC Comm_Train = %.3g\n', d_dac.comm_train);
end
if isfield(d_ac,'comm_train')
    fprintf('AC  Comm_Train = %.3g\n', d_ac.comm_train);
end
fprintf('\nInterpretation:\n');
fprintf('  DAC tracking error: x_i(k) - mean_i u_i(k). In this GP code, Pi is fixed, so the reference is constant.\n');
fprintf('  DAC disagreement  : x_i(k) - mean_i x_i(k). This checks whether agents agree.\n');
fprintf('  AC average error  : x_i(k) - mean_i x_i(0). The initial state is the local statistic Pi, not zero.\n');
