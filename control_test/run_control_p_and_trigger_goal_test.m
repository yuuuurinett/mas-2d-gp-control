function run_control_p_and_trigger_goal_test()
% Reproducible 30 s control-task test for the advisor's target:
% decreasing absolute P changes and decreasing per-agent/per-point DAC events.

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
Cleanup = onCleanup(@() restore_environment(EnvironmentNames,OldValues)); %#ok<NASGU>

setenv('CONTROL_SIM_END_TIME','30');
setenv('CONTROL_IP_NUM_INDUCING_POINTS','400');
setenv('CONTROL_DAC_TRIGGER_EPSILON','0.005');
setenv('CONTROL_DAC_TRACKING_TOL','0.005');
setenv('CONTROL_DAC_CONSENSUS_TOL','0.005');
setenv('CONTROL_DAC_TERMINATION_CRITERION','tracking');
setenv('CONTROL_DAC_KAPPA_P','1');

ResultBaseName = 'control_task_p_and_trigger_current30';
run_simulation_inducing_point( ...
    'poe',ResultFolder,ResultBaseName,true,0.2,0.1,0,42);

ResultPath = fullfile(ResultFolder,[ResultBaseName '.mat']);
Result = load(ResultPath,'t_set','AbsolutePChangeRMS', ...
    'DACBroadcastCountPerSecond','NumInducingPoints', ...
    'EarlyPChange','LatePChange');

ValidPChange = isfinite(Result.AbsolutePChangeRMS) & ...
    Result.AbsolutePChangeRMS > 0;
UpdateTimes = Result.t_set(ValidPChange);
PChange = Result.AbsolutePChangeRMS(ValidPChange);

WindowCount = numel(Result.DACBroadcastCountPerSecond);
WindowCentres = (0:WindowCount-1)+0.5;

FigureHandle = figure('Color','w','Position',[100,80,1050,700]);
Layout = tiledlayout(FigureHandle,2,1,'TileSpacing','compact', ...
    'Padding','compact');

nexttile(Layout,1);
semilogy(UpdateTimes,PChange,'o-','LineWidth',1.35,'MarkerSize',3);
grid on; box on;
xlim([0,Result.t_set(end)]);
xlabel('Time (s)');
ylabel('Absolute P change RMS');
title('Control task: absolute reference-signal change');

nexttile(Layout,2);
bar(WindowCentres,Result.DACBroadcastCountPerSecond,1);
grid on; box on;
xlim([0,Result.t_set(end)]);
xlabel('One-second window centre (s)');
ylabel('Triggers / agent / inducing point');
title('Control task: DAC trigger count in each 1 s window');

title(Layout,'Online LocalGP–IP-DAC control task', ...
    'FontSize',14,'FontWeight','bold');

FigurePath = fullfile(ResultFolder, ...
    'control_task_p_and_trigger_current30.png');
exportgraphics(FigureHandle,FigurePath,'Resolution',180, ...
    'ContentType','image');
close(FigureHandle);

ComparisonWindowCount = min(5,WindowCount);
EarlyTriggerMean = mean(Result.DACBroadcastCountPerSecond( ...
    1:ComparisonWindowCount));
LateTriggerMean = mean(Result.DACBroadcastCountPerSecond( ...
    end-ComparisonWindowCount+1:end));
TriggerReductionPercent = 100*(1-LateTriggerMean/EarlyTriggerMean);

SummaryPath = fullfile(ResultFolder, ...
    'control_task_p_and_trigger_current30_summary.txt');
FileID = fopen(SummaryPath,'w');
assert(FileID >= 0,'Could not create summary file.');
FileCleanup = onCleanup(@() fclose(FileID)); %#ok<NASGU>
fprintf(FileID,'Control task: POE/IP-DAC, 30 s\n');
fprintf(FileID,'Agents: 6\n');
fprintf(FileID,'Inducing points: %d\n',Result.NumInducingPoints);
fprintf(FileID,'Early P-change mean: %.9g\n',Result.EarlyPChange);
fprintf(FileID,'Late P-change mean: %.9g\n',Result.LatePChange);
fprintf(FileID,'Early trigger mean (first %d windows): %.9g\n', ...
    ComparisonWindowCount,EarlyTriggerMean);
fprintf(FileID,'Late trigger mean (last %d windows): %.9g\n', ...
    ComparisonWindowCount,LateTriggerMean);
fprintf(FileID,'Trigger reduction: %.3f %%\n',TriggerReductionPercent);
fprintf(FileID,'\nWindow start (s)\tTriggers / agent / inducing point\n');
for WindowNr = 1:WindowCount
    fprintf(FileID,'%d\t%.9g\n',WindowNr-1, ...
        Result.DACBroadcastCountPerSecond(WindowNr));
end

fprintf('\nControl-task goal-test summary\n');
fprintf('Early/late absolute P-change mean: %.6g -> %.6g\n', ...
    Result.EarlyPChange,Result.LatePChange);
fprintf('Early/late trigger-window mean: %.6g -> %.6g\n', ...
    EarlyTriggerMean,LateTriggerMean);
fprintf('Trigger reduction: %.2f%%\n',TriggerReductionPercent);
fprintf('Figure: %s\n',FigurePath);
fprintf('Summary: %s\n',SummaryPath);
end

function restore_environment(Names,Values)
for EnvironmentNr = 1:numel(Names)
    setenv(Names{EnvironmentNr},Values{EnvironmentNr});
end
end
