%% plot_consensus_tracking_with_disagreement_values.m
% Clean consensus diagnostics for all aggregation methods.
%
% This script does exactly three things:
%   1) Prints numerical IP-DAC inter-agent disagreement values:
%        mean max |xi_i(k) - mean_j xi_j(k)|
%      This checks whether agents agree with each other. It is NOT tracking error.
%
%   2) Plots IP-DAC representative tracking trajectories:
%        xi_i(k) tracks average reference signal mean_i p_i(k)
%
%   3) Plots IP-AC representative tracking trajectories:
%        xi_i(k) converges to initial average mean_i xi_i(0)
%
% No separate figure is generated for the DAC disagreement curve.

clc; close all;

%% settings
DatasetName = 'KIN40K';
methods = {'poe','gpoe','moe','bcm','rbcm'};
train_ratio = 0.4;
seed = 1;

tr_tag = round(train_ratio * 100);
SaveFolder = fullfile('Result','Dataset',DatasetName);

%% Load results
results = struct([]);
valid_methods = {};

for mi = 1:numel(methods)
    method = methods{mi};

    dac_file = fullfile(SaveFolder, sprintf('%s_tr%d_mc%d.mat', method, tr_tag, seed));
    ac_file  = fullfile(SaveFolder, sprintf('%s_ac_tr%d_mc%d.mat', method, tr_tag, seed));

    if ~exist(dac_file,'file')
        warning('DAC result file not found for %s: %s', upper(method), dac_file);
        continue;
    end
    if ~exist(ac_file,'file')
        warning('AC result file not found for %s: %s', upper(method), ac_file);
        continue;
    end

    d_dac = load(dac_file);
    d_ac  = load(ac_file);

    idx = numel(results) + 1;
    valid_methods{idx} = upper(method); %#ok<SAGROW>

    results(idx).method = method;
    results(idx).method_name = upper(method);
    results(idx).dac_file = dac_file;
    results(idx).ac_file = ac_file;
    results(idx).d_dac = d_dac;
    results(idx).d_ac = d_ac;

    results(idx).dac_disagreement = get_curve(d_dac, 'conv_dac_disagreement_curve', dac_file);
    results(idx).dac_tracking     = get_curve(d_dac, 'conv_dac_tracking_curve', dac_file);
    results(idx).ac_tracking      = get_curve(d_ac,  'conv_ac_avg_error_curve', ac_file);

    results(idx).dac_comm_train = get_field_or_default(d_dac, 'comm_train', NaN);
    results(idx).ac_comm_train  = get_field_or_default(d_ac,  'comm_train', NaN);

    results(idx).conv_r_dac = get_field_or_default(d_dac, 'conv_r', 1);
    results(idx).conv_m_dac = get_field_or_default(d_dac, 'conv_m', 1);
    results(idx).conv_r_ac  = get_field_or_default(d_ac,  'conv_r', 1);
    results(idx).conv_m_ac  = get_field_or_default(d_ac,  'conv_m', 1);
end

if isempty(results)
    error(['No valid result files found. Run the dataset first, for example:\n' ...
           '  run_inducingpoint_dataset(''%s'',''all'',%.2f,%d)'], ...
           DatasetName, train_ratio, seed);
end

%% 1) Numerical IP-DAC inter-agent disagreement values
% This value corresponds to:
%   mean max |xi_i(k) - mean_j xi_j(k)|

fprintf('\n====================================================================\n');
fprintf('IP-DAC inter-agent disagreement values [%s, train=%d%%, seed=%d]\n', ...
    DatasetName, tr_tag, seed);
fprintf('Quantity: mean max |xi_i(k) - mean_j xi_j(k)|\n');
fprintf('Meaning : agreement among agents; this is NOT tracking error.\n');
fprintf('====================================================================\n');
fprintf('%-8s %14s %14s %12s %14s\n', ...
    'Method', 'Initial', 'Final', 'Iterations', 'Comm_Train');
fprintf('%-8s %14s %14s %12s %14s\n', ...
    '------', '-------', '-----', '----------', '----------');

for idx = 1:numel(results)
    curve = results(idx).dac_disagreement;

    if isempty(curve)
        fprintf('%-8s %14s %14s %12s %14s\n', ...
            results(idx).method_name, 'N/A', 'N/A', 'N/A', 'N/A');
    else
        fprintf('%-8s %14.6e %14.6e %12d %14.4f\n', ...
            results(idx).method_name, curve(1), curve(end), numel(curve), results(idx).dac_comm_train);
    end
end
fprintf('====================================================================\n\n');

%% 2) Figure: DAC representative tracking trajectories
% Representative plot only: one information dimension r and one inducing point m.
figure('Name','Figure 1: IP-DAC representative tracking','Color','w','Position',[120,120,1200,760]);
tiledlayout(numel(results),1,'TileSpacing','compact','Padding','compact');

for idx = 1:numel(results)
    d_dac = results(idx).d_dac;
    nexttile;

    if isfield(d_dac,'dac_state_hist') && isfield(d_dac,'dac_ref_hist') && ~isempty(d_dac.dac_state_hist)
        plot(d_dac.dac_state_hist, 'LineWidth', 1.0); hold on;
        plot(d_dac.dac_ref_hist, 'k--', 'LineWidth', 1.4);
        grid on;
        ylabel('xi_i(k)');
        title(sprintf('%s: IP-DAC tracking of average reference signal mean_i p_i(k), r=%d, m=%d', ...
            results(idx).method_name, results(idx).conv_r_dac, results(idx).conv_m_dac), ...
            'Interpreter','none');
        if idx == 1
            legend_entries = create_agent_legend(size(d_dac.dac_state_hist,2), 'Average reference mean_i p_i(k)');
            legend(legend_entries, 'Location','bestoutside','Interpreter','none');
        end
    else
        text(0.1,0.5,sprintf('No DAC representative tracking history for %s', results(idx).method_name), ...
            'Interpreter','none');
        axis off;
    end
end
xlabel('Iteration');
sgtitle('IP-DAC representative tracking: each agent state xi_i(k) tracks average reference signal', ...
    'Interpreter','none');

%% 3) Figure: AC representative tracking trajectories
% Representative plot only: one information dimension r and one inducing point m.
figure('Name','Figure 2: IP-AC representative tracking','Color','w','Position',[140,140,1200,760]);
tiledlayout(numel(results),1,'TileSpacing','compact','Padding','compact');

for idx = 1:numel(results)
    d_ac = results(idx).d_ac;
    nexttile;

    if isfield(d_ac,'ac_state_hist') && isfield(d_ac,'ac_ref_hist') && ~isempty(d_ac.ac_state_hist)
        plot(d_ac.ac_state_hist, 'LineWidth', 1.0); hold on;
        plot(d_ac.ac_ref_hist, 'k--', 'LineWidth', 1.4);
        grid on;
        ylabel('xi_i(k)');
        title(sprintf('%s: IP-AC convergence to initial average mean_i xi_i(0), r=%d, m=%d', ...
            results(idx).method_name, results(idx).conv_r_ac, results(idx).conv_m_ac), ...
            'Interpreter','none');
        if idx == 1
            legend_entries = create_agent_legend(size(d_ac.ac_state_hist,2), 'Initial average mean_i xi_i(0)');
            legend(legend_entries, 'Location','bestoutside','Interpreter','none');
        end
    else
        text(0.1,0.5,sprintf('No AC representative tracking history for %s', results(idx).method_name), ...
            'Interpreter','none');
        axis off;
    end
end
xlabel('Iteration');
sgtitle('IP-AC representative tracking: each agent state xi_i(k) converges to initial average', ...
    'Interpreter','none');

%% Done
fprintf('Consensus tracking plots generated for %d methods: ', numel(results));
fprintf('%s ', valid_methods{:});
fprintf('\n');

%% Local helper functions
function curve = get_curve(s, field_name, file_name)
    if isfield(s, field_name)
        curve = s.(field_name);
        curve = curve(:);
    else
        curve = [];
        warning('%s not found in %s', field_name, file_name);
    end
end

function value = get_field_or_default(s, field_name, default_value)
    if isfield(s, field_name)
        value = s.(field_name);
    else
        value = default_value;
    end
end

function legend_entries = create_agent_legend(num_agents, ref_label)
    legend_entries = cell(1, num_agents + 1);
    for ii = 1:num_agents
        legend_entries{ii} = sprintf('Agent %d: xi_%d', ii, ii);
    end
    legend_entries{num_agents + 1} = ref_label;
end
