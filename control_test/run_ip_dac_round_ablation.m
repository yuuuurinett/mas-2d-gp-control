function Summary = run_ip_dac_round_ablation(do_run,round_values,m_value,simulation_end_time)
%RUN_IP_DAC_ROUND_ABLATION Sensitivity to DAC rounds per physical step.

if nargin < 1 || isempty(do_run), do_run = true; end
if nargin < 2 || isempty(round_values), round_values = [1 5 10 20 50]; end
if nargin < 3 || isempty(m_value), m_value = 600; end
if nargin < 4 || isempty(simulation_end_time), simulation_end_time = 10; end

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result','ip_dac_round_ablation');
if ~isfolder(result_folder), mkdir(result_folder); end

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_DAC_STEPS_PER_PHYSICAL_STEP','CONTROL_IP_PROJECTION_DIAGNOSTICS'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment_rounds(env_names,old_values)); %#ok<NASGU>
setenv(env_names{1},num2str(simulation_end_time,'%.15g'));
setenv(env_names{2},num2str(m_value));
setenv(env_names{4},'1');
seed = 42;

if do_run
    for idx = 1:numel(round_values)
        rounds = round_values(idx);
        setenv(env_names{3},num2str(rounds));
        fprintf('\n[DAC round ablation %d/%d] rounds=%d\n', ...
            idx,numel(round_values),rounds);
        run_simulation_inducing_point('poe',result_folder, ...
            sprintf('ip_dac_R%d_M%d',rounds,m_value),true,[],[],[],seed);
    end
end

n = numel(round_values);
RoundsPerPhysicalStep = round_values(:);
FinalTrackingError = nan(n,1);
MeanTrackingError = nan(n,1);
MeanPredictionError = nan(n,1);
MeanDirectPredictionError = nan(n,1);
MeanProjectionGap = nan(n,1);
MeanDACAverageTrackingError = nan(n,1);
FinalDACAverageTrackingError = nan(n,1);
AgentPhysicalActivationRatio = nan(n,1);
PointTriggersPerAgentPerPointPerSecond = nan(n,1);
PointTriggersPerAgent = nan(n,1);
ElapsedTimeSeconds = nan(n,1);

colors = lines(n);
fig = figure('Color','w','Position',[100 80 1150 720]);
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
track_ax = nexttile(layout); set(track_ax,'YScale','log'); hold(track_ax,'on');
pred_ax = nexttile(layout); hold(pred_ax,'on');

for idx = 1:n
    rounds = RoundsPerPhysicalStep(idx);
    data = load(fullfile(result_folder, ...
        sprintf('ip_dac_R%d_M%d.mat',rounds,m_value)));
    t = data.t_set(:);
    tracking = data.TrackingError_vector(:);
    prediction = data.prediction_error_norm_vector(:);
    plot(track_ax,t,max(tracking,eps),'Color',colors(idx,:), ...
        'LineWidth',1.6,'DisplayName',sprintf('%d rounds',rounds));
    plot(pred_ax,t,prediction,'Color',colors(idx,:), ...
        'LineWidth',1.4,'DisplayName',sprintf('%d rounds',rounds));

    FinalTrackingError(idx) = tracking(end);
    MeanTrackingError(idx) = trapz(t,tracking)/(t(end)-t(1));
    MeanPredictionError(idx) = mean(prediction,'omitnan');
    MeanDirectPredictionError(idx) = ...
        mean(data.direct_prediction_error_norm_vector,'omitnan');
    MeanProjectionGap(idx) = mean(data.ip_projection_gap_norm_vector,'omitnan');
    MeanDACAverageTrackingError(idx) = ...
        mean(data.dac_tracking_error_set,'omitnan');
    last_valid = find(isfinite(data.dac_tracking_error_set),1,'last');
    FinalDACAverageTrackingError(idx) = data.dac_tracking_error_set(last_valid);
    event_set = squeeze(any(data.dac_trigger_count_per_agent_point_set>0,2));
    AgentPhysicalActivationRatio(idx) = mean(event_set(:));
    PointTriggersPerAgent(idx) = mean(data.dac_event_count_per_agent);
    PointTriggersPerAgentPerPointPerSecond(idx) = ...
        PointTriggersPerAgent(idx)/(m_value*simulation_end_time);
    ElapsedTimeSeconds(idx) = data.elapsed_time;
end

xlabel(track_ax,'Time (s)'); ylabel(track_ax,'Tracking error');
title(track_ax,'Closed-loop tracking'); legend(track_ax,'Location','best'); grid(track_ax,'on'); box(track_ax,'on');
xlabel(pred_ax,'Time (s)'); ylabel(pred_ax,'Prediction error');
title(pred_ax,'Unknown-dynamics prediction'); legend(pred_ax,'Location','best'); grid(pred_ax,'on'); box(pred_ax,'on');

dac_ax = nexttile(layout); hold(dac_ax,'on');
plot(dac_ax,RoundsPerPhysicalStep,MeanDACAverageTrackingError,'LineWidth',2, ...
    'DisplayName','Mean');
plot(dac_ax,RoundsPerPhysicalStep,FinalDACAverageTrackingError,'LineWidth',2, ...
    'DisplayName','Final');
xlabel(dac_ax,'DAC rounds / physical step'); ylabel(dac_ax,'DAC average-tracking error');
legend(dac_ax,'Location','best'); grid(dac_ax,'on'); box(dac_ax,'on');

trigger_ax = nexttile(layout); yyaxis(trigger_ax,'left');
plot(trigger_ax,RoundsPerPhysicalStep,PointTriggersPerAgentPerPointPerSecond, ...
    'LineWidth',2); ylabel(trigger_ax,'Triggers / agent / point / s');
yyaxis(trigger_ax,'right');
plot(trigger_ax,RoundsPerPhysicalStep,AgentPhysicalActivationRatio,'LineWidth',2);
ylabel(trigger_ax,'Agent physical activation ratio'); ylim(trigger_ax,[0 1.05]);
xlabel(trigger_ax,'DAC rounds / physical step'); grid(trigger_ax,'on'); box(trigger_ax,'on');
title(layout,sprintf('IP-ET-DAC round sensitivity, M = %d',m_value));
exportgraphics(fig,fullfile(result_folder,'round_ablation.png'),'Resolution',220);
savefig(fig,fullfile(result_folder,'round_ablation.fig'));

Summary = table(RoundsPerPhysicalStep,FinalTrackingError,MeanTrackingError, ...
    MeanPredictionError,MeanDirectPredictionError,MeanProjectionGap, ...
    MeanDACAverageTrackingError,FinalDACAverageTrackingError, ...
    AgentPhysicalActivationRatio,PointTriggersPerAgent, ...
    PointTriggersPerAgentPerPointPerSecond,ElapsedTimeSeconds);
writetable(Summary,fullfile(result_folder,'summary.csv'));
save(fullfile(result_folder,'summary.mat'),'Summary','round_values', ...
    'm_value','simulation_end_time','seed');
disp(Summary);
end

function restore_environment_rounds(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
