%% 60 s control-task diagnostic with six agents kept separate

clear; clc;

ProjectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(ProjectRoot));
ResultFolder = fullfile(ProjectRoot,'control_test','result','Diagnostics');
if ~isfolder(ResultFolder)
    mkdir(ResultFolder);
end

EnvironmentNames = { ...
    'CONTROL_SIM_END_TIME', ...
    'CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_DAC_TRIGGER_EPSILON', ...
    'CONTROL_DAC_TRACKING_TOL', ...
    'CONTROL_DAC_CONSENSUS_TOL', ...
    'CONTROL_DAC_TERMINATION_CRITERION', ...
    'CONTROL_DAC_KAPPA_P'};
OldValues = cellfun(@getenv,EnvironmentNames,'UniformOutput',false);
EnvironmentCleanup = onCleanup( ...
    @() restore_environment(EnvironmentNames,OldValues)); %#ok<NASGU>

setenv('CONTROL_SIM_END_TIME','60');
setenv('CONTROL_IP_NUM_INDUCING_POINTS','400');
setenv('CONTROL_DAC_TRIGGER_EPSILON','0.05');
setenv('CONTROL_DAC_TRACKING_TOL','0.85');
setenv('CONTROL_DAC_CONSENSUS_TOL','0.005');
setenv('CONTROL_DAC_TERMINATION_CRITERION','tracking');
setenv('CONTROL_DAC_KAPPA_P','5');

ResultBaseName = 'control_task_agent_diagnostics_tau085_kappa5_eps005_60';
ResultPath = fullfile(ResultFolder,[ResultBaseName '.mat']);
if ~isfile(ResultPath)
    run_simulation_inducing_point( ...
        'poe',ResultFolder,ResultBaseName,true,0.2,0.1,0,42);
end

Result = load(ResultPath,'t_set','projection_update_set', ...
    'projection_step_change_per_agent_set', ...
    'dac_tracking_error_per_agent_set','dac_broadcast_count_set', ...
    'NumInducingPoints','DACTriggerEpsilon','DACTrackingTolerance','Kappa_P');

AgentQuantity = size(Result.dac_broadcast_count_set,1);
AgentColors = lines(AgentQuantity);
PhysicalTime = Result.t_set(1:end-1);
OnlineUpdateMask = Result.projection_update_set(1:end-1) > 0;
OnlineUpdateTime = PhysicalTime(OnlineUpdateMask);
PChange = Result.projection_step_change_per_agent_set(:,OnlineUpdateMask);
TrackingError = Result.dac_tracking_error_per_agent_set(:,OnlineUpdateMask);
TriggerCount = Result.dac_broadcast_count_set(:,OnlineUpdateMask);
SimulationEndTime = Result.t_set(end);

%% Six-agent P change
PChangeFigure = figure('Color','w','Position',[80,80,1250,540]);
hold on;
for AgentNr = 1:AgentQuantity
    semilogy(OnlineUpdateTime,PChange(AgentNr,:),'-', ...
        'Color',AgentColors(AgentNr,:),'LineWidth',0.9, ...
        'DisplayName',sprintf('Agent %d',AgentNr));
end
grid on; box on; xlim([0,SimulationEndTime]);
xlabel('Time (s)'); ylabel('RMS of P_i(k)-P_i(k-1)');
title('Control task: absolute P change for six agents');
legend('Location','eastoutside');
PChangePath = fullfile(ResultFolder,[ResultBaseName '_p_change_6_agents.png']);
exportgraphics(PChangeFigure,PChangePath,'Resolution',180);
close(PChangeFigure);

%% Six-agent final DAC tracking error after every online update
TrackingFigure = figure('Color','w','Position',[80,80,1250,540]);
hold on;
for AgentNr = 1:AgentQuantity
    semilogy(OnlineUpdateTime,TrackingError(AgentNr,:),'-', ...
        'Color',AgentColors(AgentNr,:),'LineWidth',0.9, ...
        'DisplayName',sprintf('Agent %d',AgentNr));
end
yline(Result.DACTrackingTolerance,'k--','Tracking tolerance 0.85', ...
    'LabelHorizontalAlignment','left','HandleVisibility','off');
grid on; box on; xlim([0,SimulationEndTime]);
xlabel('Time (s)');
ylabel('max_m ||Xi_{i,m}-mean_j(P_{j,m})||_2');
title('Control task: final DAC tracking error for six agents');
legend('Location','eastoutside');
TrackingPath = fullfile(ResultFolder, ...
    [ResultBaseName '_tracking_error_6_agents.png']);
exportgraphics(TrackingFigure,TrackingPath,'Resolution',180);
close(TrackingFigure);

%% Six-agent broadcasts at every 0.1 s update
TriggerFigure = figure('Color','w','Position',[80,80,1250,540]);
hold on;
for AgentNr = 1:AgentQuantity
    plot(OnlineUpdateTime,TriggerCount(AgentNr,:),'-', ...
        'Color',AgentColors(AgentNr,:),'LineWidth',0.8, ...
        'DisplayName',sprintf('Agent %d',AgentNr));
end
grid on; box on; xlim([0,SimulationEndTime]);
xlabel('Time (s)'); ylabel('Broadcasts per inducing point');
title('Control task: DAC broadcasts at every 0.1 s update');
legend('Location','eastoutside');
TriggerPath = fullfile(ResultFolder, ...
    [ResultBaseName '_trigger_count_0p1s.png']);
exportgraphics(TriggerFigure,TriggerPath,'Resolution',180);
close(TriggerFigure);

%% Simple-MAS-style one-second star plot
SecondEdges = 0:1:SimulationEndTime;
BroadcastPerSecond = zeros(AgentQuantity,numel(SecondEdges)-1);
for SecondNr = 1:numel(SecondEdges)-1
    WindowMask = OnlineUpdateTime >= SecondEdges(SecondNr) & ...
        OnlineUpdateTime < SecondEdges(SecondNr+1);
    BroadcastPerSecond(:,SecondNr) = sum(TriggerCount(:,WindowMask),2);
end

StarFigure = figure('Color','w','Position',[80,80,1250,520]);
hold on;
for AgentNr = 1:AgentQuantity
    StarQuantity = round(BroadcastPerSecond(AgentNr,:));
    for SecondNr = 1:numel(StarQuantity)
        if StarQuantity(SecondNr) > 0
            StarTime = linspace(SecondEdges(SecondNr)+0.08, ...
                SecondEdges(SecondNr+1)-0.08,StarQuantity(SecondNr));
            plot(StarTime,AgentNr*ones(size(StarTime)),'*', ...
                'LineStyle','none','Color',AgentColors(AgentNr,:), ...
                'MarkerSize',3.5);
        end
    end
end
grid on; box on; xlim([0,SimulationEndTime]);
ylim([0.5,AgentQuantity+0.5]); yticks(1:AgentQuantity);
xlabel('Time (s)'); ylabel('Agent');
title(['Average DAC broadcasts per agent per inducing point ', ...
    '(one star = one average broadcast)']);
StarPath = fullfile(ResultFolder,[ResultBaseName '_trigger_stars.png']);
exportgraphics(StarFigure,StarPath,'Resolution',180);
close(StarFigure);

%% First/last 5 s numerical comparison without median
EarlyMask = OnlineUpdateTime < 5;
LateMask = OnlineUpdateTime >= SimulationEndTime-5;
EarlyP = mean(PChange(:,EarlyMask),2,'omitnan');
LateP = mean(PChange(:,LateMask),2,'omitnan');
EarlyTracking = mean(TrackingError(:,EarlyMask),2);
LateTracking = mean(TrackingError(:,LateMask),2);
EarlyTrigger = mean(BroadcastPerSecond(:,1:5),2);
LateTrigger = mean(BroadcastPerSecond(:,end-4:end),2);

SummaryPath = fullfile(ResultFolder,[ResultBaseName '_summary.md']);
FileID = fopen(SummaryPath,'w');
assert(FileID >= 0,'Cannot create the 60 s diagnostic summary.');
FileCleanup = onCleanup(@() fclose(FileID)); %#ok<NASGU>
fprintf(FileID,'# Control task: 60 s six-agent diagnostic\n\n');
fprintf(FileID,['Tracking tolerance = 0.85, Kappa_P = 5, epsilon = 0.05, ', ...
    '400 inducing points. All values are arithmetic means.\n\n']);
fprintf(FileID,['| Agent | P change, first 5 s | P change, last 5 s | ', ...
    'Tracking error, first 5 s | Tracking error, last 5 s | ', ...
    'Triggers, first 5 s | Triggers, last 5 s |\n']);
fprintf(FileID,'|---:|---:|---:|---:|---:|---:|---:|\n');
for AgentNr = 1:AgentQuantity
    fprintf(FileID,'| %d | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g |\n', ...
        AgentNr,EarlyP(AgentNr),LateP(AgentNr), ...
        EarlyTracking(AgentNr),LateTracking(AgentNr), ...
        EarlyTrigger(AgentNr),LateTrigger(AgentNr));
end
fprintf(FileID,'| Mean | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g |\n', ...
    mean(EarlyP),mean(LateP),mean(EarlyTracking),mean(LateTracking), ...
    mean(EarlyTrigger),mean(LateTrigger));

fprintf('\n60 s means: first 5 s -> last 5 s\n');
fprintf('P change: %.6g -> %.6g\n',mean(EarlyP),mean(LateP));
fprintf('DAC tracking error: %.6g -> %.6g\n', ...
    mean(EarlyTracking),mean(LateTracking));
fprintf('Triggers: %.6g -> %.6g\n', ...
    mean(EarlyTrigger),mean(LateTrigger));
fprintf('P plot: %s\nTracking plot: %s\nTrigger plot: %s\n', ...
    PChangePath,TrackingPath,TriggerPath);
fprintf('Star plot: %s\nMarkdown: %s\n',StarPath,SummaryPath);

function restore_environment(Names,Values)
for EnvironmentNr = 1:numel(Names)
    setenv(Names{EnvironmentNr},Values{EnvironmentNr});
end
end
