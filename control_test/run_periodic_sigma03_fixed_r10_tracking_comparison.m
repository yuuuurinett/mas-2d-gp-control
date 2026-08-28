function Summary = run_periodic_sigma03_fixed_r10_tracking_comparison(do_run)
%RUN_PERIODIC_SIGMA03_FIXED_R10_TRACKING_COMPARISON Fair 30-s PoE comparison.

if nargin < 1 || isempty(do_run), do_run = true; end
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result', ...
    'periodic_sigma03_10slots_per_0p1s_T30');
if ~isfolder(result_folder), mkdir(result_folder); end
point_file = fullfile(repo_root,'result','ip_reachable_m600', ...
    'reachable_uniform_M600.mat');

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_IP_INDUCING_POINT_FILE','CONTROL_IP_PROJECTION_DIAGNOSTICS', ...
    'CONTROL_DAC_STEPS_PER_PHYSICAL_STEP','CONTROL_DAC_TRIGGER_EPSILON', ...
    'CONTROL_AC_ITERATION_POLICY','CONTROL_AC_FIXED_ITERATIONS', ...
    'CONTROL_AC_BROADCAST_TRIGGER','CONTROL_AC_PERIODIC_SIGMA', ...
    'CONTROL_AC_TRIGGER_DIAGNOSTICS','CONTROL_TP_DAC_STEPS_PER_PHYSICAL_STEP', ...
    'CONTROL_TP_DAC_INITIAL_ROUNDS','CONTROL_TP_QUERY_UPDATE_INTERVAL', ...
    'CONTROL_TP_AC_ITERATION_POLICY','CONTROL_TP_AC_FIXED_ROUNDS'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment_local(env_names,old_values)); %#ok<NASGU>

setenv('CONTROL_SIM_END_TIME','30');
setenv('CONTROL_IP_NUM_INDUCING_POINTS','600');
setenv('CONTROL_IP_INDUCING_POINT_FILE',point_file);
setenv('CONTROL_IP_PROJECTION_DIAGNOSTICS','0');
% DAC evolves once per 0.01-s physical step, hence 10 updates in each
% 0.1-s online-learning interval. AC performs the same 10 updates as an
% inner batch whenever P is refreshed.
setenv('CONTROL_DAC_STEPS_PER_PHYSICAL_STEP','1');
setenv('CONTROL_DAC_TRIGGER_EPSILON','0.005');
setenv('CONTROL_AC_ITERATION_POLICY','fixed');
setenv('CONTROL_AC_FIXED_ITERATIONS','10');
setenv('CONTROL_AC_BROADCAST_TRIGGER','petc_cache');
setenv('CONTROL_AC_PERIODIC_SIGMA','0.3');
setenv('CONTROL_AC_TRIGGER_DIAGNOSTICS','0');
setenv('CONTROL_TP_DAC_STEPS_PER_PHYSICAL_STEP','1');
setenv('CONTROL_TP_DAC_INITIAL_ROUNDS','1');
setenv('CONTROL_TP_QUERY_UPDATE_INTERVAL','0.1');
setenv('CONTROL_TP_AC_ITERATION_POLICY','fixed');
setenv('CONTROL_TP_AC_FIXED_ROUNDS','10');

paths.ip_dac = fullfile(result_folder,'ip_dac_poe_10slots.mat');
paths.ip_ac = fullfile(result_folder,'ip_ac_petc_sigma03_R10.mat');
paths.tp_dac = fullfile(result_folder,'tp_dac_poe_10slots.mat');
paths.tp_ac = fullfile(result_folder,'tp_ac_poe_R10.mat');
paths.cen = fullfile(result_folder,'cen.mat');
paths.nbr = fullfile(result_folder,'nbr.mat');
paths.local = fullfile(result_folder,'local_offline_peragent350_seeded.mat');
paths.exact = fullfile(result_folder,'exact.mat');

if do_run
    run_if_missing(paths.ip_dac,@() run_simulation_inducing_point( ...
        'poe',result_folder,'ip_dac_poe_10slots',false,[],[],[],42));
    run_if_missing(paths.ip_ac,@() run_simulation_inducing_point( ...
        'poe_ac',result_folder,'ip_ac_petc_sigma03_R10',false,[],[],[],42));
    run_if_missing(paths.tp_dac,@() run_simulation_test_point( ...
        'poe',result_folder,'tp_dac_poe_10slots',false,42,0));
    run_if_missing(paths.tp_ac,@() run_simulation_test_point( ...
        'poe_ac',result_folder,'tp_ac_poe_R10',false,42,0));
    run_if_missing(paths.cen,@() run_simulation_cen( ...
        'poe',result_folder,'cen',false,42,0));
    run_if_missing(paths.nbr,@() run_simulation_nbr( ...
        'poe',result_folder,'nbr',false,42,0));
    % Offline Local uses 350 fixed samples per agent, matching shared code.
    run_if_missing(paths.local,@() run_simulation_test_point( ...
        'local_offline',result_folder,'local_offline_peragent350_seeded', ...
        false,42,350));
    run_if_missing(paths.exact,@() run_simulation_test_point( ...
        'exact',result_folder,'exact',false,42,0));
end

specs = [ ...
    struct('key','ip_dac','label','IP-DAC (10 updates / 0.1 s)','comm','ip_dac'); ...
    struct('key','ip_ac','label','IP-AC (PETC, sigma=0.3, R=10)','comm','ip_ac'); ...
    struct('key','tp_dac','label','TP-DAC (10 updates / 0.1 s)','comm','tp'); ...
    struct('key','tp_ac','label','TP-AC (R=10)','comm','tp'); ...
    struct('key','cen','label','CEN','comm','none'); ...
    struct('key','nbr','label','NBR','comm','none'); ...
    struct('key','local','label','Local (offline, 350 / agent)','comm','none'); ...
    struct('key','exact','label','Exact','comm','none')];

n = numel(specs);
Method = strings(n,1);
MaxTrackingErrorAfter10s = nan(n,1);
BroadcastsPerAgent = nan(n,1);
history = cell(n,1);
time = cell(n,1);
for idx = 1:n
    d = load(paths.(specs(idx).key));
    Method(idx) = specs(idx).label;
    time{idx} = d.t_set(:);
    if isfield(d,'vartheta_all_set')
        history{idx} = sqrt(sum(d.vartheta_all_set.^2,1)).';
    else
        history{idx} = d.TrackingError_vector(:);
    end
    stable_mask = time{idx} >= 10;
    MaxTrackingErrorAfter10s(idx) = max( ...
        history{idx}(stable_mask),[],'omitnan');
    switch specs(idx).comm
        case 'ip_dac'
            BroadcastsPerAgent(idx) = d.dac_broadcasts_per_agent;
        case 'ip_ac'
            BroadcastsPerAgent(idx) = d.ac_broadcasts_per_agent;
        case 'tp'
            BroadcastsPerAgent(idx) = mean(d.tp_broadcast_count_per_agent);
    end
end

fig = figure('Visible','off','Color','w','Position',[80 50 1100 800]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax_full = nexttile(tl); hold(ax_full,'on'); grid(ax_full,'on'); box(ax_full,'on');
ax_zoom = nexttile(tl); hold(ax_zoom,'on'); grid(ax_zoom,'on'); box(ax_zoom,'on');
colors = [0.55 0.20 0.62; 0.84 0.20 0.15; 0.18 0.50 0.25; ...
          0.90 0.50 0.10; 0.25 0.55 0.55; 0.55 0.38 0.20; ...
          0.45 0.45 0.45; 0.05 0.05 0.05];
styles = {'-','-','--','--','-.','-.',':','-'};
for idx = 1:n
    plot(ax_full,time{idx},history{idx},'Color',colors(idx,:), ...
        'LineStyle',styles{idx},'LineWidth',1.35, ...
        'DisplayName',char(Method(idx)));
    plot(ax_zoom,time{idx},history{idx},'Color',colors(idx,:), ...
        'LineStyle',styles{idx},'LineWidth',1.35, ...
        'HandleVisibility','off');
end
xlim(ax_full,[0 30]); ylabel(ax_full,'Tracking error');
title(ax_full,'Full trajectory');
legend(ax_full,'Location','eastoutside');
xline(ax_full,10,'k:','HandleVisibility','off');
xlim(ax_zoom,[10 30]);
stable_plot_max = max(MaxTrackingErrorAfter10s,[],'omitnan');
ylim(ax_zoom,[0.05,1.05*stable_plot_max]);
xlabel(ax_zoom,'Time (s)'); ylabel(ax_zoom,'Tracking error');
title(ax_zoom,'Evaluation window: t >= 10 s');
title(tl,'PoE tracking comparison: 10 consensus updates per 0.1 s');
output_png = fullfile(result_folder,'tracking_error_fixed_R10.png');
print(fig,output_png,'-dpng','-r200');
savefig(fig,strrep(output_png,'.png','.fig'));
close(fig);

Summary = table(Method,MaxTrackingErrorAfter10s,BroadcastsPerAgent);
writetable(Summary,fullfile(result_folder,'tracking_summary.csv'));
disp(Summary);
end

function run_if_missing(path,runner)
if ~isfile(path), runner(); end
end

function restore_environment_local(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
