%% main_BatchTest_FINAL_PAPER
% Final paper-style batch runner + summary exporter
% Tables:
%   1) Final tracking error of all online GP baselines
%   2) Online learning / communication statistics
%   3) Markdown export for direct report use

clc; close all;

%% Configuration
CurrentMode    = 'all';       % 'poe'|'gpoe'|'moe'|'bcm'|'rbcm' or 'all'
TestType       = 'all';       % 'test'|'inducing'|'cen'|'nbr'|'all'
use_formation  = true;        % true = formation, false = no formation

% Set these false when you only want to reload existing .mat files and export tables.
force_no_run = true;

% File suffix
form_tag = 'formation';
if ~use_formation
    form_tag = 'noformation';
end

%% Paths
SaveFolder_Test     = fullfile('Result', 'Test_Point');
SaveFolder_Inducing = fullfile('Result', 'Inducing_Point');
SaveFolder_CEN      = fullfile('Result', 'CEN');
SaveFolder_NBR      = fullfile('Result', 'NBR');
SaveFolder_Figures  = fullfile('Result', 'Figures');
SaveFolder_Tables   = fullfile('Result', 'Tables');

for f = {SaveFolder_Test, SaveFolder_Inducing, SaveFolder_CEN, SaveFolder_NBR, SaveFolder_Figures, SaveFolder_Tables}
    if ~exist(f{1}, 'dir')
        mkdir(f{1});
    end
end

%% Mode lists
Modes_dac      = {'poe','gpoe','moe','bcm','rbcm'};
Modes_ac       = {'poe_ac','gpoe_ac','moe_ac','bcm_ac','rbcm_ac'};
Modes_baseline = {'local','exact'};

if strcmpi(CurrentMode, 'all')
    AllModes_dac = Modes_dac;
else
    AllModes_dac = {lower(CurrentMode)};
end

%% 1. Run simulations
run_inducing = ismember(lower(TestType), {'inducing','all'});
run_test     = ismember(lower(TestType), {'test','all'});
run_cen      = ismember(lower(TestType), {'cen','all'});
run_nbr      = ismember(lower(TestType), {'nbr','all'});


%{
if force_no_run
    run_inducing = false;
    run_test     = false;
    run_cen      = false;
    run_nbr      = false;
end
%}

if run_inducing
    AllModes_ind = [Modes_dac, Modes_ac, Modes_baseline];
    fprintf('\n======== Inducing-point aggregation: IP-DAC + IP-AC [%s] ========\n', form_tag);
    for m = 1:numel(AllModes_ind)
        fname = sprintf('%s_%s', AllModes_ind{m}, form_tag);
        fprintf('[%d/%d] %s\n', m, numel(AllModes_ind), fname);
        run_simulation_inducing_point(AllModes_ind{m}, SaveFolder_Inducing, fname, use_formation);
    end
end

if run_test
    AllModes_test = [Modes_dac, Modes_ac, Modes_baseline];
    fprintf('\n======== Test-point aggregation: TP-DAC + TP-AC [%s] ========\n', form_tag);
    for m = 1:numel(AllModes_test)
        fname = sprintf('%s_%s', AllModes_test{m}, form_tag);
        fprintf('[%d/%d] %s\n', m, numel(AllModes_test), fname);
        run_simulation_test_point(AllModes_test{m}, SaveFolder_Test, fname, use_formation);
    end
end

if run_cen
    fprintf('\n======== Centralized aggregation: CEN [%s] ========\n', form_tag);
    for m = 1:numel(AllModes_dac)
        fname = sprintf('cen_%s_%s', AllModes_dac{m}, form_tag);
        fprintf('[%d/%d] %s\n', m, numel(AllModes_dac), fname);
        run_simulation_cen(AllModes_dac{m}, SaveFolder_CEN, fname, use_formation);
    end
end

if run_nbr
    fprintf('\n======== Neighbor aggregation: NBR [%s] ========\n', form_tag);
    for m = 1:numel(AllModes_dac)
        fname = sprintf('nbr_%s_%s', AllModes_dac{m}, form_tag);
        fprintf('[%d/%d] %s\n', m, numel(AllModes_dac), fname);
        run_simulation_nbr(AllModes_dac{m}, SaveFolder_NBR, fname, use_formation);
    end
end

%% 2. Load time axis
ref_candidates = {
    fullfile(SaveFolder_Inducing, sprintf('poe_%s.mat', form_tag));
    fullfile(SaveFolder_Test,     sprintf('poe_%s.mat', form_tag));
    fullfile(SaveFolder_CEN,      sprintf('cen_poe_%s.mat', form_tag));
    fullfile(SaveFolder_NBR,      sprintf('nbr_poe_%s.mat', form_tag))};

t_set = [];
for i = 1:numel(ref_candidates)
    if exist(ref_candidates{i}, 'file')
        tmp = load(ref_candidates{i}, 't_set');
        if isfield(tmp, 't_set')
            t_set = tmp.t_set;
            break;
        end
    end
end

if isempty(t_set)
    fprintf('No result file found. Nothing to summarize.\n');
    return;
end

N = numel(t_set);
load_err = @(folder, fname) load_tracking_error(folder, fname, N);

%% 3. Load tracking errors
Err_IP_DAC = nan(numel(Modes_dac), N);
Err_IP_AC  = nan(numel(Modes_ac),  N);
Err_TP_DAC = nan(numel(Modes_dac), N);
Err_TP_AC  = nan(numel(Modes_ac),  N);
Err_CEN    = nan(numel(Modes_dac), N);
Err_NBR    = nan(numel(Modes_dac), N);

for m = 1:numel(Modes_dac)
    Err_IP_DAC(m,:) = load_err(SaveFolder_Inducing, sprintf('%s_%s', Modes_dac{m}, form_tag));
    Err_TP_DAC(m,:) = load_err(SaveFolder_Test,     sprintf('%s_%s', Modes_dac{m}, form_tag));
    Err_CEN(m,:)    = load_err(SaveFolder_CEN,      sprintf('cen_%s_%s', Modes_dac{m}, form_tag));
    Err_NBR(m,:)    = load_err(SaveFolder_NBR,      sprintf('nbr_%s_%s', Modes_dac{m}, form_tag));
end

for m = 1:numel(Modes_ac)
    Err_IP_AC(m,:) = load_err(SaveFolder_Inducing, sprintf('%s_%s', Modes_ac{m}, form_tag));
    Err_TP_AC(m,:) = load_err(SaveFolder_Test,     sprintf('%s_%s', Modes_ac{m}, form_tag));
end

Err_Local = load_err(SaveFolder_Test, sprintf('local_%s', form_tag));
Err_Exact = load_err(SaveFolder_Test, sprintf('exact_%s', form_tag));

%% 4. Load online-learning and communication statistics
Stats = struct();

for m = 1:numel(Modes_dac)
    method = Modes_dac{m};
    acmethod = Modes_ac{m};

    f_ip_dac = fullfile(SaveFolder_Inducing, sprintf('%s_%s.mat', method, form_tag));
    f_ip_ac  = fullfile(SaveFolder_Inducing, sprintf('%s_%s.mat', acmethod, form_tag));
    f_tp_dac = fullfile(SaveFolder_Test,     sprintf('%s_%s.mat', method, form_tag));
    f_tp_ac  = fullfile(SaveFolder_Test,     sprintf('%s_%s.mat', acmethod, form_tag));
    f_cen    = fullfile(SaveFolder_CEN,      sprintf('cen_%s_%s.mat', method, form_tag));
    f_nbr    = fullfile(SaveFolder_NBR,      sprintf('nbr_%s_%s.mat', method, form_tag));

    Stats(m).method = method; 
    Stats(m).IP_DAC = load_gp_stats(f_ip_dac);
    Stats(m).IP_AC  = load_gp_stats(f_ip_ac);
    Stats(m).TP_DAC = load_gp_stats(f_tp_dac);
    Stats(m).TP_AC  = load_gp_stats(f_tp_ac);
    Stats(m).CEN    = load_gp_stats(f_cen);
    Stats(m).NBR    = load_gp_stats(f_nbr);
end

%% 5. Print final tracking-error table
fprintf('\n%s\n', repmat('=', 1, 104));
fprintf('  Final tracking error ||e(T)||  [%s]\n', form_tag);
fprintf('  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s\n', ...
    'Method','IP-DAC','IP-AC','TP-DAC','TP-AC','CEN','NBR','Local','Exact');
fprintf('  %s\n', repmat('-', 1, 100));

for m = 1:numel(Modes_dac)
    fprintf('  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s\n', ...
        upper(Modes_dac{m}), ...
        fmt_num(Err_IP_DAC(m,end), '%.4f'), ...
        fmt_num(Err_IP_AC(m,end),  '%.4f'), ...
        fmt_num(Err_TP_DAC(m,end), '%.4f'), ...
        fmt_num(Err_TP_AC(m,end),  '%.4f'), ...
        fmt_num(Err_CEN(m,end),    '%.4f'), ...
        fmt_num(Err_NBR(m,end),    '%.4f'), ...
        fmt_num(Err_Local(end),    '%.4f'), ...
        fmt_num(Err_Exact(end),    '%.4f'));
end
fprintf('%s\n', repmat('=', 1, 104));

%% 6. Print inducing-point consensus communication table
fprintf('\n%s\n', repmat('=', 1, 118));
fprintf('  IP consensus broadcast triggers: average broadcasts / agent  [%s]\n', form_tag);
fprintf('  %-8s  %-18s  %-18s\n', 'Method','IP-DAC','IP-AC');
fprintf('  %s\n', repmat('-', 1, 114));

for m = 1:numel(Modes_dac)
    fprintf('  %-8s  %-18s  %-18s\n', ...
        upper(Modes_dac{m}), ...
        fmt_num(Stats(m).IP_DAC.dac_comm_ag, '%.2f'), ...
        fmt_num(Stats(m).IP_AC.ac_comm_ag,   '%.2f'));
end
fprintf('%s\n', repmat('=', 1, 118));

%% 7. Export Markdown tables
md_file = fullfile(SaveFolder_Tables, sprintf('online_learning_summary_%s.md', form_tag));
export_markdown_summary(md_file, form_tag, Modes_dac, Err_IP_DAC, Err_IP_AC, Err_TP_DAC, Err_TP_AC, Err_CEN, Err_NBR, Err_Local, Err_Exact, Stats);
fprintf('\nMarkdown summary exported to:\n%s\n', md_file);

%% 8. Plot tracking curves
lw = 1.5;
for m = 1:numel(Modes_dac)
    ac_name = [Modes_dac{m}, '_ac'];
    ac_idx  = find(strcmp(Modes_ac, ac_name));

    fig = figure('Name', upper(Modes_dac{m}), 'Color', 'w');
    hold on; grid on; box on;
    set(gca, 'YScale', 'log', 'FontSize', 11, 'FontName', 'Times New Roman');

    if ~all(isnan(Err_IP_DAC(m,:)))
        plot(t_set, Err_IP_DAC(m,:), 'b-',  'LineWidth', lw, 'DisplayName', 'IP-DAC');
    end
    if ~isempty(ac_idx) && ~all(isnan(Err_IP_AC(ac_idx,:)))
        plot(t_set, Err_IP_AC(ac_idx,:), 'b--', 'LineWidth', lw, 'DisplayName', 'IP-AC');
    end
    if ~all(isnan(Err_TP_DAC(m,:)))
        plot(t_set, Err_TP_DAC(m,:), 'r-',  'LineWidth', lw, 'DisplayName', 'TP-DAC');
    end
    if ~isempty(ac_idx) && ~all(isnan(Err_TP_AC(ac_idx,:)))
        plot(t_set, Err_TP_AC(ac_idx,:), 'r--', 'LineWidth', lw, 'DisplayName', 'TP-AC');
    end
    if ~all(isnan(Err_CEN(m,:)))
        plot(t_set, Err_CEN(m,:), 'g-',  'LineWidth', lw, 'DisplayName', 'CEN');
    end
    if ~all(isnan(Err_NBR(m,:)))
        plot(t_set, Err_NBR(m,:), 'm-',  'LineWidth', lw, 'DisplayName', 'NBR');
    end
    if ~all(isnan(Err_Local))
        plot(t_set, Err_Local, 'k--', 'LineWidth', 1, 'DisplayName', 'Local');
    end
    if ~all(isnan(Err_Exact))
        plot(t_set, Err_Exact, 'k-',  'LineWidth', 1, 'DisplayName', 'Exact');
    end

    ylabel('$\|e\|$', 'Interpreter', 'latex', 'FontSize', 13);
    xlabel('$t$ (s)', 'Interpreter', 'latex', 'FontSize', 13);
    title(sprintf('Tracking Error - %s [%s]', upper(Modes_dac{m}), form_tag), ...
        'FontSize', 12, 'FontName', 'Times New Roman');
    legend('Location', 'northeast', 'FontSize', 9, 'NumColumns', 2);
    xlim([t_set(1), t_set(end)]);

    fname_fig = sprintf('Comparison_%s_%s', upper(Modes_dac{m}), form_tag);
    saveas(fig, fullfile(SaveFolder_Figures, [fname_fig, '.png']));
    savefig(fig, fullfile(SaveFolder_Figures, [fname_fig, '.fig']));

    dac_result_file = fullfile(SaveFolder_Inducing, ...
        sprintf('%s_%s.mat', Modes_dac{m}, form_tag));
    ac_result_file = fullfile(SaveFolder_Inducing, ...
        sprintf('%s_%s.mat', ac_name, form_tag));
    plot_agent_tracking_errors(dac_result_file, ac_result_file, ...
        Modes_dac{m}, form_tag, SaveFolder_Figures);
    plot_consensus_broadcasts(dac_result_file, ac_result_file, ...
        Modes_dac{m}, form_tag, SaveFolder_Figures);
end

fprintf('\nDone. Figures saved to %s\n', SaveFolder_Figures);

%% Local helper functions
function err = load_tracking_error(folder, fname, N)
    fpath = fullfile(folder, [fname, '.mat']);
    err = nan(1, N);

    if exist(fpath, 'file')
        d = load(fpath, 'TrackingError_vector');
        if isfield(d, 'TrackingError_vector')
            len = min(numel(d.TrackingError_vector), N);
            err(1:len) = d.TrackingError_vector(1:len);
        end
    end
end

function st = load_gp_stats(fpath)
    st = struct();
    st.online_avg   = NaN;
    st.dac_comm_ag  = NaN;
    st.ac_comm_ag   = NaN;

    if ~exist(fpath, 'file')
        return;
    end

    d = load(fpath);

    % Fixed time-triggered online-learning updates.
    if isfield(d, 'online_trigger_count')
        st.online_avg = mean(d.online_trigger_count(:));
    elseif isfield(d, 'online_update_count')
        st.online_avg = mean(d.online_update_count(:));
    end

    % DAC communication: average agent-level broadcasts.
    if isfield(d, 'dac_broadcasts_per_agent')
        st.dac_comm_ag = d.dac_broadcasts_per_agent;
    elseif isfield(d, 'dac_trigger_count_per_agent_point')
        st.dac_comm_ag = d.dac_trigger_count_per_agent_point;
    elseif isfield(d, 'trigger_count_per_agent_point')
        st.dac_comm_ag = d.trigger_count_per_agent_point;
    elseif isfield(d, 'dac_total_trigger_count') && isfield(d, 'NumInducingPoints')
        st.dac_comm_ag = mean(d.dac_total_trigger_count(:)) / d.NumInducingPoints;
    elseif isfield(d, 'dac_total_trigger_count') && isfield(d, 'InducingPoints_Coordinates')
        point_count = size(d.InducingPoints_Coordinates,2);
        st.dac_comm_ag = mean(d.dac_total_trigger_count(:)) / point_count;
    elseif isfield(d, 'total_trigger_count') && isfield(d, 'NumInducingPoints')
        st.dac_comm_ag = mean(d.total_trigger_count(:)) / d.NumInducingPoints;
    end

    % AC communication: average agent-level broadcasts.
    if isfield(d, 'ac_broadcasts_per_agent')
        st.ac_comm_ag = d.ac_broadcasts_per_agent;
    elseif isfield(d, 'ac_trigger_count_per_agent_point')
        st.ac_comm_ag = d.ac_trigger_count_per_agent_point;
    elseif isfield(d, 'ac_trigger_count_per_agent')
        st.ac_comm_ag = d.ac_trigger_count_per_agent;
    elseif isfield(d, 'ac_total_trigger_count')
        st.ac_comm_ag = mean(d.ac_total_trigger_count(:));
    end
end

function plot_agent_tracking_errors(dac_file, ac_file, method, form_tag, output_folder)
    dac = load_agent_tracking_history(dac_file);
    ac = load_agent_tracking_history(ac_file);
    if isempty(dac.error) && isempty(ac.error)
        return;
    end

    agent_count = max(size(dac.error,1), size(ac.error,1));
    fig = figure('Name', sprintf('%s agent tracking errors', upper(method)), ...
        'Color', 'w');
    layout = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', ...
        'Padding', 'compact');
    for agent_i = 1:agent_count
        ax = nexttile(layout);
        hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
        set(ax, 'YScale', 'log', 'FontName', 'Times New Roman');
        if agent_i <= size(dac.error,1)
            plot(ax, dac.time, dac.error(agent_i,:), 'b-', ...
                'LineWidth', 1.2, 'DisplayName', 'IP-DAC');
        end
        if agent_i <= size(ac.error,1)
            plot(ax, ac.time, ac.error(agent_i,:), 'r--', ...
                'LineWidth', 1.2, 'DisplayName', 'IP-AC');
        end
        title(ax, sprintf('Agent %d', agent_i));
        xlabel(ax, 't (s)'); ylabel(ax, '$\|\vartheta_i(t)\|$', ...
            'Interpreter', 'latex');
        if agent_i == 1
            legend(ax, 'Location', 'best');
        end
    end
    title(layout, sprintf('Tracking error of each agent - %s [%s]', ...
        upper(method), form_tag));
    base_name = sprintf('AgentTracking_%s_%s', upper(method), form_tag);
    saveas(fig, fullfile(output_folder, [base_name, '.png']));
    savefig(fig, fullfile(output_folder, [base_name, '.fig']));
end

function history = load_agent_tracking_history(fpath)
    history = struct('error', [], 'time', []);
    if ~isfile(fpath)
        return;
    end
    data = load(fpath, 'vartheta_all_set', 't_set');
    if ~isfield(data, 'vartheta_all_set') || ~isfield(data, 't_set')
        return;
    end
    state_dimension = 4;
    agent_count = size(data.vartheta_all_set,1) / state_dimension;
    if agent_count ~= floor(agent_count)
        return;
    end
    reshaped_error = reshape(data.vartheta_all_set, ...
        state_dimension, agent_count, []);
    history.error = squeeze(vecnorm(reshaped_error,2,1));
    history.time = data.t_set(1:size(history.error,2));
end

function plot_consensus_broadcasts(dac_file, ac_file, method, form_tag, output_folder)
    dac = load_broadcast_history(dac_file, 'dac_broadcast_count_set');
    ac = load_broadcast_history(ac_file, 'ac_broadcast_count_set');
    if isempty(dac.counts) && isempty(ac.counts)
        return;
    end

    fig = figure('Name', sprintf('%s consensus broadcasts', upper(method)), ...
        'Color', 'w');
    layout = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', ...
        'Padding', 'compact');

    ax1 = nexttile(layout);
    hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');
    plot_broadcast_instances(ax1, dac, 0, [0 0.4470 0.7410], 'DAC');
    agent_offset = 0;
    if ~isempty(dac.counts)
        agent_offset = size(dac.counts,1) + 1;
    end
    plot_broadcast_instances(ax1, ac, agent_offset, [0.8500 0.3250 0.0980], 'AC');
    xlabel(ax1, 't (s)'); ylabel(ax1, 'Agent');
    title(ax1, 'Broadcast trigger instances');

    ax2 = nexttile(layout);
    hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
    if ~isempty(dac.counts)
        plot(ax2, dac.time, mean(cumsum(dac.counts,2),1), 'b-', ...
            'LineWidth', 1.5, 'DisplayName', 'IP-DAC');
    end
    if ~isempty(ac.counts)
        plot(ax2, ac.time, mean(cumsum(ac.counts,2),1), 'r--', ...
            'LineWidth', 1.5, 'DisplayName', 'IP-AC');
    end
    xlabel(ax2, 't (s)'); ylabel(ax2, 'Average cumulative broadcasts / agent');
    title(ax2, 'Communication growth over simulation time');
    legend(ax2, 'Location', 'northwest');

    title(layout, sprintf('%s consensus communication [%s]', ...
        upper(method), form_tag));
    base_name = sprintf('ConsensusTriggers_%s_%s', upper(method), form_tag);
    saveas(fig, fullfile(output_folder, [base_name, '.png']));
    savefig(fig, fullfile(output_folder, [base_name, '.fig']));
end

function history = load_broadcast_history(fpath, field_name)
    history = struct('counts', [], 'time', []);
    if ~isfile(fpath)
        return;
    end
    data = load(fpath, field_name, 't_set');
    if ~isfield(data, field_name) || ~isfield(data, 't_set')
        return;
    end
    history.counts = data.(field_name);
    sample_count = size(history.counts,2);
    history.time = data.t_set(1:sample_count);
end

function plot_broadcast_instances(ax, history, agent_offset, color, label_prefix)
    if isempty(history.counts)
        return;
    end
    agent_count = size(history.counts,1);
    for agent_i = 1:agent_count
        event_idx = find(history.counts(agent_i,:) > 0);
        if isempty(event_idx)
            continue;
        end
        event_sizes = 12 + 8*log1p(history.counts(agent_i,event_idx));
        scatter(ax, history.time(event_idx), ...
            (agent_offset+agent_i)*ones(size(event_idx)), event_sizes, ...
            color, 'x', 'DisplayName', sprintf('%s agent %d', label_prefix, agent_i));
    end
end

function s = fmt_num(x, fmt)
    if isempty(x) || isnan(x)
        s = '-';
    else
        s = sprintf(fmt, x);
    end
end

function export_markdown_summary(md_file, form_tag, Modes_dac, Err_IP_DAC, Err_IP_AC, Err_TP_DAC, Err_TP_AC, Err_CEN, Err_NBR, Err_Local, Err_Exact, Stats)
    fid = fopen(md_file, 'w');
    if fid == -1
        error('Cannot open markdown file for writing: %s', md_file);
    end

    fprintf(fid, '# Online learning comparison under %s control\n\n', form_tag);
    fprintf(fid, 'All GP-based baselines use online data-triggered local GP updates. In formation control, neighboring GP models are evaluated at the controlled agent state `x_i`, not at the neighbor state `x_j`.\n\n');

    fprintf(fid, '## Final tracking error\n\n');
    fprintf(fid, '| Method | IP-DAC | IP-AC | TP-DAC | TP-AC | CEN | NBR | Local | Exact |\n');
    fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|---:|---:|\n');
    for m = 1:numel(Modes_dac)
        fprintf(fid, '| %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', ...
            upper(Modes_dac{m}), ...
            fmt_num(Err_IP_DAC(m,end), '%.4f'), ...
            fmt_num(Err_IP_AC(m,end),  '%.4f'), ...
            fmt_num(Err_TP_DAC(m,end), '%.4f'), ...
            fmt_num(Err_TP_AC(m,end),  '%.4f'), ...
            fmt_num(Err_CEN(m,end),    '%.4f'), ...
            fmt_num(Err_NBR(m,end),    '%.4f'), ...
            fmt_num(Err_Local(end),    '%.4f'), ...
            fmt_num(Err_Exact(end),    '%.4f'));
    end

    fprintf(fid, '\n## IP consensus communication triggers\n\n');
    fprintf(fid, '| Method | IP-DAC broadcasts / agent | IP-AC broadcasts / agent |\n');
    fprintf(fid, '|---|---:|---:|\n');
    for m = 1:numel(Modes_dac)
        fprintf(fid, '| %s | %s | %s |\n', ...
            upper(Modes_dac{m}), ...
            fmt_num(Stats(m).IP_DAC.dac_comm_ag, '%.2f'), ...
            fmt_num(Stats(m).IP_AC.ac_comm_ag,   '%.2f'));
    end

    fprintf(fid, '\n## Notes\n\n');
    fprintf(fid, '- `IP` denotes inducing-point aggregation.\n');
    fprintf(fid, '- `TP` denotes test-point aggregation.\n');
    fprintf(fid, '- `DAC` uses dynamic average consensus.\n');
    fprintf(fid, '- `AC` uses static average consensus.\n');
    fprintf(fid, '- `CEN` uses centralized aggregation.\n');
    fprintf(fid, '- `NBR` uses neighbor-only aggregation.\n');
    fprintf(fid, '- `Exact` is a non-GP reference and therefore has no online GP update.\n');

    fclose(fid);
end
