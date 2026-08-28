function Summary = run_ip_m_ablation(do_run,m_values,simulation_end_time,dac_rounds)
%RUN_IP_M_ABLATION Fixed-seed inducing-point count sensitivity study.
%
% Only M changes. Inducing-point generation, controller, online sampling,
% IP-DAC rounds, and event-trigger parameters remain identical.

if nargin < 1 || isempty(do_run), do_run = true; end
if nargin < 2 || isempty(m_values), m_values = [100 200 400 600 800]; end
if nargin < 3 || isempty(simulation_end_time), simulation_end_time = 10; end
if nargin < 4 || isempty(dac_rounds), dac_rounds = 10; end

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result','ip_m_ablation');
if ~isfolder(result_folder), mkdir(result_folder); end

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_DAC_STEPS_PER_PHYSICAL_STEP', ...
    'CONTROL_IP_PROJECTION_DIAGNOSTICS'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment_m(env_names,old_values)); %#ok<NASGU>
setenv(env_names{1},num2str(simulation_end_time,'%.15g'));
setenv(env_names{3},num2str(dac_rounds));
setenv(env_names{4},'1');

seed = 42;
if do_run
    for m_idx = 1:numel(m_values)
        M = m_values(m_idx);
        setenv(env_names{2},num2str(M));
        fprintf('\n[M ablation %d/%d] M=%d\n',m_idx,numel(m_values),M);
        run_simulation_inducing_point('poe',result_folder, ...
            sprintf('ip_et_dac_poe_M%d',M),true,[],[],[],seed);
    end
end

method_count = numel(m_values);
M = m_values(:);
FinalTrackingError = nan(method_count,1);
FullTimeMeanTrackingError = nan(method_count,1);
MeanPredictionError = nan(method_count,1);
MeanDirectPredictionError = nan(method_count,1);
MeanProjectionGap = nan(method_count,1);
PointTriggersPerAgent = nan(method_count,1);
PointTriggersPerAgentPerPoint = nan(method_count,1);
PointTriggersPerAgentPerPointPerSecond = nan(method_count,1);
ElapsedTimeSeconds = nan(method_count,1);

colors = lines(method_count);
tracking_fig = figure('Color','w','Position',[100 100 1100 650]);
tracking_ax = axes(tracking_fig);
set(tracking_ax,'YScale','log'); hold(tracking_ax,'on');

for m_idx = 1:method_count
    m_i = M(m_idx);
    path_i = fullfile(result_folder,sprintf('ip_et_dac_poe_M%d.mat',m_i));
    if ~isfile(path_i), error('Missing M-ablation result: %s',path_i); end
    data = load(path_i);
    t = data.t_set(:);
    e = data.TrackingError_vector(:);
    semilogy(tracking_ax,t,max(e,eps),'LineWidth',1.8, ...
        'Color',colors(m_idx,:),'DisplayName',sprintf('M = %d',m_i));

    FinalTrackingError(m_idx) = e(end);
    FullTimeMeanTrackingError(m_idx) = trapz(t,e)/(t(end)-t(1));
    MeanPredictionError(m_idx) = ...
        mean(data.prediction_error_norm_vector,'omitnan');
    MeanDirectPredictionError(m_idx) = ...
        mean(data.direct_prediction_error_norm_vector,'omitnan');
    MeanProjectionGap(m_idx) = ...
        mean(data.ip_projection_gap_norm_vector,'omitnan');
    PointTriggersPerAgent(m_idx) = mean(data.dac_event_count_per_agent);
    PointTriggersPerAgentPerPoint(m_idx) = ...
        PointTriggersPerAgent(m_idx)/m_i;
    PointTriggersPerAgentPerPointPerSecond(m_idx) = ...
        PointTriggersPerAgentPerPoint(m_idx)/simulation_end_time;
    ElapsedTimeSeconds(m_idx) = data.elapsed_time;
end

xlabel(tracking_ax,'Time (s)');
ylabel(tracking_ax,'Overall tracking error ||\vartheta(t)||_2');
title(tracking_ax,'IP-ET-DAC sensitivity to inducing-point count');
legend(tracking_ax,'Location','best'); grid(tracking_ax,'on'); box(tracking_ax,'on');
xlim(tracking_ax,[0 simulation_end_time]);
exportgraphics(tracking_fig,fullfile(result_folder,'tracking_error_vs_M.png'), ...
    'Resolution',220);
savefig(tracking_fig,fullfile(result_folder,'tracking_error_vs_M.fig'));

metric_fig = figure('Color','w','Position',[120 120 1100 720]);
tiledlayout(metric_fig,2,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(M,FullTimeMeanTrackingError,'LineWidth',2);
xlabel('M'); ylabel('Mean tracking error'); grid on; box on;
nexttile; hold on;
plot(M,MeanDirectPredictionError,'LineWidth',2,'DisplayName','Direct PoE');
plot(M,MeanPredictionError,'LineWidth',2,'DisplayName','IP reconstruction');
xlabel('M'); ylabel('Mean prediction error'); legend('Location','best'); grid on; box on;
nexttile; plot(M,MeanProjectionGap,'LineWidth',2);
xlabel('M'); ylabel('Mean projection gap'); grid on; box on;
nexttile; plot(M,PointTriggersPerAgentPerPointPerSecond,'LineWidth',2);
xlabel('M'); ylabel('Triggers / agent / point / s'); grid on; box on;
exportgraphics(metric_fig,fullfile(result_folder,'metrics_vs_M.png'), ...
    'Resolution',220);
savefig(metric_fig,fullfile(result_folder,'metrics_vs_M.fig'));

Summary = table(M,FinalTrackingError,FullTimeMeanTrackingError, ...
    MeanPredictionError,MeanDirectPredictionError,MeanProjectionGap, ...
    PointTriggersPerAgent,PointTriggersPerAgentPerPoint, ...
    PointTriggersPerAgentPerPointPerSecond,ElapsedTimeSeconds);
writetable(Summary,fullfile(result_folder,'summary.csv'));
save(fullfile(result_folder,'summary.mat'),'Summary','m_values', ...
    'simulation_end_time','dac_rounds','seed');
disp(Summary);
end

function restore_environment_m(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
