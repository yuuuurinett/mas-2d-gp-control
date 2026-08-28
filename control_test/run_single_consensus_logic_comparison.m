function Summary = run_single_consensus_logic_comparison( ...
    do_run,simulation_end_time,num_inducing_points,dac_rounds)
%RUN_SINGLE_CONSENSUS_LOGIC_COMPARISON Fixed-seed PoE AC/DAC validation.
%
% Physical sampling: 0.01 s. Online LocalGP update: every 0.1 s.
% Each agent evaluates only GP_i(x_i). TP-AC restarts from the resulting
% p-by-N snapshot; TP-DAC preserves one p-by-N Zeta state across snapshots.
% IP-ET-AC restarts only after each inducing projection refresh.
% IP-ET-DAC preserves Zeta and broadcast snapshots for the full run.

if nargin < 1 || isempty(do_run), do_run = true; end
if nargin < 2 || isempty(simulation_end_time), simulation_end_time = 10; end
if nargin < 3 || isempty(num_inducing_points), num_inducing_points = 400; end
if nargin < 4 || isempty(dac_rounds), dac_rounds = 10; end

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result','single_consensus_logic');
if ~isfolder(result_folder), mkdir(result_folder); end

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_DAC_STEPS_PER_PHYSICAL_STEP', ...
    'CONTROL_TP_DAC_STEPS_PER_PHYSICAL_STEP'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment(env_names,old_values)); %#ok<NASGU>
setenv(env_names{1},num2str(simulation_end_time,'%.15g'));
setenv(env_names{2},num2str(num_inducing_points));
setenv(env_names{3},num2str(dac_rounds));
setenv(env_names{4},num2str(dac_rounds));

seed = 42;
if do_run
    run_simulation_test_point('exact',result_folder,'exact',true,seed,0);
    run_simulation_test_point('local',result_folder,'local',true,seed,0);
    offline_total = 350;
    offline_allocation = floor(offline_total/6)*ones(6,1);
    offline_allocation(1:mod(offline_total,6)) = ...
        offline_allocation(1:mod(offline_total,6))+1;
    run_simulation_nbr('poe_offline',result_folder,'offline_nbr_poe', ...
        true,seed,offline_allocation);
    run_simulation_cen('poe',result_folder,'cen_poe',true,seed,0);
    run_simulation_test_point('poe_ac',result_folder,'tp_ac_poe',true,seed,0);
    run_simulation_test_point('poe',result_folder,'tp_dac_poe',true,seed,0);
    run_simulation_inducing_point('poe_ac',result_folder,'ip_et_ac_poe', ...
        true,[],[],[],seed);
    run_simulation_inducing_point('poe',result_folder,'ip_et_dac_poe', ...
        true,[],[],[],seed);
end

specs = [ ...
    struct('key','exact','label','Exact'); ...
    struct('key','local','label','Local'); ...
    struct('key','offline_nbr_poe','label','Offline-NBR-PoE (350 total)'); ...
    struct('key','cen_poe','label','CEN-PoE'); ...
    struct('key','tp_ac_poe','label','TP-AC-PoE'); ...
    struct('key','tp_dac_poe','label',sprintf('TP-DAC-PoE (%d rounds)',dac_rounds)); ...
    struct('key','ip_et_ac_poe','label','IP-ET-AC-PoE'); ...
    struct('key','ip_et_dac_poe','label',sprintf('IP-ET-DAC-PoE (%d rounds)',dac_rounds))];

method_count = numel(specs);
Method = strings(method_count,1);
FinalTrackingError = nan(method_count,1);
FullTimeMeanTrackingError = nan(method_count,1);
ElapsedTimeSeconds = nan(method_count,1);
OnlineSamplesPerAgent = nan(method_count,1);
InitialOfflineSamplesPerAgent = nan(method_count,1);
CommunicationEventsPerAgent = nan(method_count,1);

colors = [0 0 0; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; ...
    0.6350 0.0780 0.1840; 0 0.4470 0.7410; ...
    0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560; ...
    0.4660 0.6740 0.1880];
styles = {'-','-','--','-','--',':','-.','-'};
fig = figure('Color','w','Position',[100 100 1100 650]);
ax = axes(fig);
set(ax,'YScale','log');
hold(ax,'on');

initial_errors = nan(method_count,1);
for method_idx = 1:method_count
    path_i = fullfile(result_folder,[specs(method_idx).key '.mat']);
    if ~isfile(path_i), error('Missing result: %s',path_i); end
    data = load(path_i);
    t = data.t_set(:);
    e = data.TrackingError_vector(:);
    if any(~isfinite(e)), error('%s is numerically invalid.',specs(method_idx).label); end

    semilogy(ax,t,max(e,eps),'LineStyle',styles{method_idx}, ...
        'Color',colors(method_idx,:),'LineWidth',1.8, ...
        'DisplayName',specs(method_idx).label);
    Method(method_idx) = specs(method_idx).label;
    initial_errors(method_idx) = e(1);
    FinalTrackingError(method_idx) = e(end);
    FullTimeMeanTrackingError(method_idx) = trapz(t,e)/(t(end)-t(1));
    if isfield(data,'elapsed_time'), ElapsedTimeSeconds(method_idx) = data.elapsed_time; end
    if isfield(data,'online_trigger_count')
        OnlineSamplesPerAgent(method_idx) = mean(data.online_trigger_count);
    elseif isfield(data,'online_update_count')
        OnlineSamplesPerAgent(method_idx) = mean(data.online_update_count);
    end
    if isfield(data,'OfflineDataQuantity')
        InitialOfflineSamplesPerAgent(method_idx) = data.OfflineDataQuantity;
    elseif isfield(data,'OfflineDataQuantity_set')
        InitialOfflineSamplesPerAgent(method_idx) = ...
            mean(data.OfflineDataQuantity_set);
    end
    if isfield(data,'tp_broadcast_count_per_agent')
        CommunicationEventsPerAgent(method_idx) = ...
            mean(data.tp_broadcast_count_per_agent);
    elseif contains(specs(method_idx).key,'ip_et_ac') && ...
            isfield(data,'ac_event_count_per_agent')
        CommunicationEventsPerAgent(method_idx) = mean(data.ac_event_count_per_agent);
    elseif contains(specs(method_idx).key,'ip_et_dac') && ...
            isfield(data,'dac_event_count_per_agent')
        CommunicationEventsPerAgent(method_idx) = mean(data.dac_event_count_per_agent);
    end
end

assert(max(initial_errors)-min(initial_errors)<1e-12, ...
    'Methods do not share the same initial condition.');
xlabel(ax,'Time (s)');
ylabel(ax,'Overall tracking error ||\vartheta(t)||_2');
title(ax,'Single-run PoE consensus comparison');
legend(ax,'Location','best'); grid(ax,'on'); box(ax,'on');
xlim(ax,[0 simulation_end_time]);
exportgraphics(fig,fullfile(result_folder,'tracking_error_comparison.png'), ...
    'Resolution',220);
savefig(fig,fullfile(result_folder,'tracking_error_comparison.fig'));

Summary = table(Method,FinalTrackingError,FullTimeMeanTrackingError, ...
    InitialOfflineSamplesPerAgent,OnlineSamplesPerAgent, ...
    CommunicationEventsPerAgent,ElapsedTimeSeconds);
writetable(Summary,fullfile(result_folder,'summary.csv'));
save(fullfile(result_folder,'summary.mat'),'Summary','simulation_end_time', ...
    'num_inducing_points','dac_rounds','seed');
disp(Summary);
end

function restore_environment(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
