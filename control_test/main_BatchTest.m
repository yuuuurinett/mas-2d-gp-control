%% main_BatchTest_FINAL_PAPER
% Final paper-style batch runner + summary exporter
% Tables:
%   1) Final tracking error of all online GP baselines
%   2) Online learning / communication statistics
%   3) Markdown export for direct report use

clc; close all;

%% Configuration
CurrentMode    = 'all';       % 'poe'|'gpoe'|'moe'|'bcm'|'rbcm' or 'all'
TestType       = 'inducing';       % 'test'|'inducing'|'cen'|'nbr'|'all'
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
%run_inducing = ismember(lower(TestType), {'inducing','all'});
%run_test     = ismember(lower(TestType), {'test','all'});
%run_cen      = ismember(lower(TestType), {'cen','all'});
%run_nbr      = ismember(lower(TestType), {'nbr','all'});


if force_no_run
    run_inducing = false;
    run_test     = false;
    run_cen      = false;
    run_nbr      = false;
end

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

%% 6. Print online-learning table
fprintf('\n%s\n', repmat('=', 1, 118));
fprintf('  Online learning ET: average triggers / agent  [%s]\n', form_tag);
fprintf('  %-8s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n', ...
    'Method','IP-DAC','IP-AC','TP-DAC','TP-AC','CEN','NBR');
fprintf('  %s\n', repmat('-', 1, 114));

for m = 1:numel(Modes_dac)
    fprintf('  %-8s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n', ...
        upper(Modes_dac{m}), ...
        fmt_num(Stats(m).IP_DAC.online_avg, '%.2f'), ...
        fmt_num(Stats(m).IP_AC.online_avg,  '%.2f'), ...
        fmt_num(Stats(m).TP_DAC.online_avg, '%.2f'), ...
        fmt_num(Stats(m).TP_AC.online_avg,  '%.2f'), ...
        fmt_num(Stats(m).CEN.online_avg,    '%.2f'), ...
        fmt_num(Stats(m).NBR.online_avg,    '%.2f'));
end
fprintf('%s\n', repmat('=', 1, 118));

%% 7. Print communication table
fprintf('\n%s\n', repmat('=', 1, 118));
fprintf('  Communication ET statistics  [%s]\n', form_tag);
fprintf('  %-8s  %-14s  %-14s  %-14s  %-14s\n', ...
    'Method','IP-DAC /pt','IP-AC /agent','TP-DAC /pt','TP-AC /agent');
fprintf('  %s\n', repmat('-', 1, 114));

for m = 1:numel(Modes_dac)
    fprintf('  %-8s  %-14s  %-14s  %-14s  %-14s\n', ...
        upper(Modes_dac{m}), ...
        fmt_num(Stats(m).IP_DAC.dac_comm_pt, '%.4f'), ...
        fmt_num(Stats(m).IP_AC.ac_comm_ag,   '%.2f'), ...
        fmt_num(Stats(m).TP_DAC.dac_comm_pt, '%.4f'), ...
        fmt_num(Stats(m).TP_AC.ac_comm_ag,   '%.2f'));
end
fprintf('%s\n', repmat('=', 1, 118));

%% 8. Export Markdown tables
md_file = fullfile(SaveFolder_Tables, sprintf('online_learning_summary_%s.md', form_tag));
export_markdown_summary(md_file, form_tag, Modes_dac, Err_IP_DAC, Err_IP_AC, Err_TP_DAC, Err_TP_AC, Err_CEN, Err_NBR, Err_Local, Err_Exact, Stats);
fprintf('\nMarkdown summary exported to:\n%s\n', md_file);

%% 9. Plot tracking curves
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
    st.dac_comm_pt  = NaN;
    st.ac_comm_ag   = NaN;

    if ~exist(fpath, 'file')
        return;
    end

    d = load(fpath);

    % Online learning data-trigger statistics
    if isfield(d, 'online_trigger_count')
        st.online_avg = mean(d.online_trigger_count(:));
    end

    % DAC communication: average triggers / agent / point
    if isfield(d, 'dac_trigger_count_per_agent_point')
        st.dac_comm_pt = d.dac_trigger_count_per_agent_point;
    elseif isfield(d, 'trigger_count_per_agent_point')
        st.dac_comm_pt = d.trigger_count_per_agent_point;
    elseif isfield(d, 'dac_total_trigger_count') && isfield(d, 'NumInducingPoints')
        st.dac_comm_pt = mean(d.dac_total_trigger_count(:)) / d.NumInducingPoints;
    elseif isfield(d, 'total_trigger_count') && isfield(d, 'NumInducingPoints')
        st.dac_comm_pt = mean(d.total_trigger_count(:)) / d.NumInducingPoints;
    end

    % AC communication: average triggers / agent
    if isfield(d, 'ac_trigger_count_per_agent')
        st.ac_comm_ag = d.ac_trigger_count_per_agent;
    elseif isfield(d, 'ac_total_trigger_count')
        st.ac_comm_ag = mean(d.ac_total_trigger_count(:));
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

    fprintf(fid, '\n## Online learning data-trigger statistics\n\n');
    fprintf(fid, '| Method | IP-DAC ET / agent | IP-AC ET / agent | TP-DAC ET / agent | TP-AC ET / agent | CEN ET / agent | NBR ET / agent |\n');
    fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|\n');
    for m = 1:numel(Modes_dac)
        fprintf(fid, '| %s | %s | %s | %s | %s | %s | %s |\n', ...
            upper(Modes_dac{m}), ...
            fmt_num(Stats(m).IP_DAC.online_avg, '%.2f'), ...
            fmt_num(Stats(m).IP_AC.online_avg,  '%.2f'), ...
            fmt_num(Stats(m).TP_DAC.online_avg, '%.2f'), ...
            fmt_num(Stats(m).TP_AC.online_avg,  '%.2f'), ...
            fmt_num(Stats(m).CEN.online_avg,    '%.2f'), ...
            fmt_num(Stats(m).NBR.online_avg,    '%.2f'));
    end

    fprintf(fid, '\n## Communication statistics\n\n');
    fprintf(fid, '| Method | IP-DAC comm. / agent / point | IP-AC comm. / agent | TP-DAC comm. / agent / point | TP-AC comm. / agent |\n');
    fprintf(fid, '|---|---:|---:|---:|---:|\n');
    for m = 1:numel(Modes_dac)
        fprintf(fid, '| %s | %s | %s | %s | %s |\n', ...
            upper(Modes_dac{m}), ...
            fmt_num(Stats(m).IP_DAC.dac_comm_pt, '%.4f'), ...
            fmt_num(Stats(m).IP_AC.ac_comm_ag,   '%.2f'), ...
            fmt_num(Stats(m).TP_DAC.dac_comm_pt, '%.4f'), ...
            fmt_num(Stats(m).TP_AC.ac_comm_ag,   '%.2f'));
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
