function Summary = run_full_mode_fixed_ac_comparison(do_run)
%RUN_FULL_MODE_FIXED_AC_COMPARISON One-seed, 30-s, no-formation comparison.

if nargin < 1 || isempty(do_run), do_run = true; end
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result','full_mode_tp_snapshot_T30');
if ~isfolder(result_folder), mkdir(result_folder); end

simulation_end_time = 30;
seed = 42;
use_formation = false;
old_end_time = getenv('CONTROL_SIM_END_TIME');
old_tp_rounds = getenv('CONTROL_TP_DAC_STEPS_PER_PHYSICAL_STEP');
old_tp_initial = getenv('CONTROL_TP_DAC_INITIAL_ROUNDS');
old_tp_query = getenv('CONTROL_TP_QUERY_UPDATE_INTERVAL');
old_tp_ac_policy = getenv('CONTROL_TP_AC_ITERATION_POLICY');
old_tp_ac_rounds = getenv('CONTROL_TP_AC_FIXED_ROUNDS');
cleanup_env = onCleanup(@() restore_environment(old_end_time,old_tp_rounds, ...
    old_tp_initial,old_tp_query,old_tp_ac_policy,old_tp_ac_rounds)); %#ok<NASGU>
setenv('CONTROL_SIM_END_TIME','30');
setenv('CONTROL_TP_DAC_STEPS_PER_PHYSICAL_STEP','1');
setenv('CONTROL_TP_DAC_INITIAL_ROUNDS','1');
setenv('CONTROL_TP_QUERY_UPDATE_INTERVAL','0.1');
setenv('CONTROL_TP_AC_ITERATION_POLICY','fixed');

baseline_specs = [ ...
    struct('key','tp_dac','label','TP-DAC','runner','tp','mode','poe','ac_rounds',0); ...
    struct('key','tp_ac10','label','TP-AC-10','runner','tp','mode','poe_ac','ac_rounds',10); ...
    struct('key','tp_ac20','label','TP-AC-20','runner','tp','mode','poe_ac','ac_rounds',20); ...
    struct('key','cen','label','CEN','runner','cen','mode','poe','ac_rounds',0); ...
    struct('key','nbr','label','NBR','runner','nbr','mode','poe','ac_rounds',0); ...
    struct('key','local','label','Local','runner','tp','mode','local','ac_rounds',0); ...
    struct('key','exact','label','Exact','runner','tp','mode','exact','ac_rounds',0)];

if do_run
    for idx = 1:numel(baseline_specs)
        spec = baseline_specs(idx);
        target = fullfile(result_folder,[spec.key '.mat']);
        if isfile(target), continue; end
        fprintf('\n[%d/%d] Running %s\n',idx,numel(baseline_specs),spec.label);
        switch spec.runner
            case 'tp'
                if spec.ac_rounds > 0
                    setenv('CONTROL_TP_AC_FIXED_ROUNDS', ...
                        num2str(spec.ac_rounds));
                end
                run_simulation_test_point(spec.mode,result_folder,spec.key, ...
                    use_formation,seed,0);
            case 'cen'
                run_simulation_cen(spec.mode,result_folder,spec.key, ...
                    use_formation,seed,0);
            case 'nbr'
                run_simulation_nbr(spec.mode,result_folder,spec.key, ...
                    use_formation,seed,0);
        end
    end
end

ip_folder = fullfile(repo_root,'result','ip_reachable_m600');
specs = [ ...
    struct('label','IP-DAC','file',fullfile(ip_folder, ...
        'poe_M600_reachable_R1_eps_0p2_kappa_1_T30.mat'),'comm','dac'); ...
    struct('label','IP-AC-10','file',fullfile(ip_folder, ...
        'poe_ac_M600_reachable_fixed_R10_T30.mat'),'comm','ac'); ...
    struct('label','IP-AC-20','file',fullfile(ip_folder, ...
        'poe_ac_M600_reachable_fixed_R20_T30.mat'),'comm','ac'); ...
    struct('label','TP-DAC','file',fullfile(result_folder,'tp_dac.mat'),'comm','periodic-dac'); ...
    struct('label','TP-AC-10','file',fullfile(result_folder,'tp_ac10.mat'),'comm','periodic-ac'); ...
    struct('label','TP-AC-20','file',fullfile(result_folder,'tp_ac20.mat'),'comm','periodic-ac'); ...
    struct('label','CEN','file',fullfile(result_folder,'cen.mat'),'comm','direct-global'); ...
    struct('label','NBR','file',fullfile(result_folder,'nbr.mat'),'comm','direct-neighbor'); ...
    struct('label','Local','file',fullfile(result_folder,'local.mat'),'comm','none'); ...
    struct('label','Exact','file',fullfile(result_folder,'exact.mat'),'comm','none')];

n = numel(specs);
Method = strings(n,1);
MeanTrackingError = nan(n,1);
FinalTrackingError = nan(n,1);
MeanPredictionError = nan(n,1);
BroadcastsPerAgent = nan(n,1);
CommunicationDefinition = strings(n,1);

fig = figure('Color','w','Position',[80 80 1100 650]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
ax.Toolbar.Visible = 'off';
colors = lines(n);
styles = {'-','--',':','-','--',':','-','-','-.','--'};
for idx = 1:n
    if ~isfile(specs(idx).file), error('Missing result: %s',specs(idx).file); end
    data = load(specs(idx).file);
    t = data.t_set(:);
    tracking = data.TrackingError_vector(:);
    Method(idx) = specs(idx).label;
    MeanTrackingError(idx) = trapz(t,tracking)/max(t(end)-t(1),eps);
    FinalTrackingError(idx) = tracking(end);
    if isfield(data,'prediction_error_norm_vector')
        MeanPredictionError(idx) = mean(data.prediction_error_norm_vector,'omitnan');
    end
    switch specs(idx).comm
        case 'dac'
            BroadcastsPerAgent(idx) = data.dac_broadcasts_per_agent;
            CommunicationDefinition(idx) = "agent-level ET broadcasts";
        case 'ac'
            BroadcastsPerAgent(idx) = data.ac_broadcasts_per_agent;
            CommunicationDefinition(idx) = "agent-level ET broadcasts";
        case 'direct-global'
            CommunicationDefinition(idx) = "ideal centralized access; not ET-counted";
        case 'direct-neighbor'
            CommunicationDefinition(idx) = "direct one-hop model access; not ET-counted";
        case 'periodic-dac'
            BroadcastsPerAgent(idx) = mean(data.tp_broadcast_count_per_agent);
            CommunicationDefinition(idx) = "periodic packed DAC broadcasts";
        case 'periodic-ac'
            BroadcastsPerAgent(idx) = mean(data.tp_broadcast_count_per_agent);
            CommunicationDefinition(idx) = "periodic packed AC broadcasts";
        otherwise
            CommunicationDefinition(idx) = "no consensus broadcast";
    end
    semilogy(ax,t,max(tracking,eps),'LineWidth',1.55, ...
        'LineStyle',styles{idx},'Color',colors(idx,:), ...
        'DisplayName',specs(idx).label);
end
xlabel(ax,'Time (s)'); ylabel(ax,'Overall tracking error ||\vartheta(t)||_2');
title(ax,'PoE control comparison: fixed-round AC and online baselines');
legend(ax,'Location','eastoutside'); xlim(ax,[0 simulation_end_time]);
set(ax,'YScale','log');
exportgraphics(fig,fullfile(result_folder,'tracking_error_all_modes.png'), ...
    'Resolution',200);

Summary = table(Method,MeanTrackingError,FinalTrackingError, ...
    MeanPredictionError,BroadcastsPerAgent,CommunicationDefinition);
writetable(Summary,fullfile(result_folder,'all_mode_summary.csv'));
disp(Summary);
end

function restore_environment(end_time,tp_rounds,tp_initial,tp_query, ...
    tp_ac_policy,tp_ac_rounds)
setenv('CONTROL_SIM_END_TIME',end_time);
setenv('CONTROL_TP_DAC_STEPS_PER_PHYSICAL_STEP',tp_rounds);
setenv('CONTROL_TP_DAC_INITIAL_ROUNDS',tp_initial);
setenv('CONTROL_TP_QUERY_UPDATE_INTERVAL',tp_query);
setenv('CONTROL_TP_AC_ITERATION_POLICY',tp_ac_policy);
setenv('CONTROL_TP_AC_FIXED_ROUNDS',tp_ac_rounds);
end
