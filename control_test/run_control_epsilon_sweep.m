%% Step 3: tune only the DAC broadcast epsilon in equation (17)
% Fixed: tracking criterion, tracking tolerance = 0.85, Kappa_P = 5.

clear; clc;

ProjectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(ProjectRoot));

ResultFolder = fullfile(ProjectRoot,'control_test','result','Diagnostics');
if ~isfolder(ResultFolder)
    mkdir(ResultFolder);
end

EpsilonSet = [0.005,0.05,0.1,0.3,0.6];
ResultFileSet = { ...
    'control_task_dac_graft_tracking_tau085_kappa5_30.mat'; ...
    'control_task_eps005_tau085_kappa5_30.mat'; ...
    'control_task_eps01_tau085_kappa5_30.mat'; ...
    'control_task_eps03_tau085_kappa5_30.mat'; ...
    'control_task_eps06_tau085_kappa5_30.mat'};

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

setenv('CONTROL_SIM_END_TIME','30');
setenv('CONTROL_IP_NUM_INDUCING_POINTS','400');
setenv('CONTROL_DAC_TRACKING_TOL','0.85');
setenv('CONTROL_DAC_CONSENSUS_TOL','0.005');
setenv('CONTROL_DAC_TERMINATION_CRITERION','tracking');
setenv('CONTROL_DAC_KAPPA_P','5');

for EpsilonNr = 1:numel(EpsilonSet)
    ResultPath = fullfile(ResultFolder,ResultFileSet{EpsilonNr});
    if isfile(ResultPath)
        fprintf('Reuse epsilon=%g: %s\n',EpsilonSet(EpsilonNr),ResultPath);
        continue;
    end

    setenv('CONTROL_DAC_TRIGGER_EPSILON',num2str(EpsilonSet(EpsilonNr),16));
    [~,ResultBaseName] = fileparts(ResultFileSet{EpsilonNr});
    fprintf('Run epsilon=%g ...\n',EpsilonSet(EpsilonNr));
    run_simulation_inducing_point( ...
        'poe',ResultFolder,ResultBaseName,true,0.2,0.1,0,42);
end

FinalControlTrackingError = zeros(size(EpsilonSet));
TotalBroadcasts = zeros(size(EpsilonSet));
EarlyTriggerMean = zeros(size(EpsilonSet));
LateTriggerMean = zeros(size(EpsilonSet));
ZeroIterationUpdates = zeros(size(EpsilonSet));
MeanDACIterations = zeros(size(EpsilonSet));

for EpsilonNr = 1:numel(EpsilonSet)
    ResultPath = fullfile(ResultFolder,ResultFileSet{EpsilonNr});
    Data = load(ResultPath,'TrackingError_vector', ...
        'dac_event_count_per_agent','DACBroadcastCountPerSecond', ...
        'dac_iteration_count_set','consensus_communication_update_set');

    FinalControlTrackingError(EpsilonNr) = Data.TrackingError_vector(end);
    TotalBroadcasts(EpsilonNr) = sum(Data.dac_event_count_per_agent);
    EarlyTriggerMean(EpsilonNr) = ...
        mean(Data.DACBroadcastCountPerSecond(1:5));
    LateTriggerMean(EpsilonNr) = ...
        mean(Data.DACBroadcastCountPerSecond(end-4:end));

    OnlineUpdateMask = ...
        Data.consensus_communication_update_set(1:end-1) > 0;
    UpdateIterationCount = ...
        Data.dac_iteration_count_set(OnlineUpdateMask);
    ZeroIterationUpdates(EpsilonNr) = sum(UpdateIterationCount == 0);
    MeanDACIterations(EpsilonNr) = mean(UpdateIterationCount);
end

SummaryPath = fullfile(ResultFolder,'control_epsilon_sweep_summary.md');
FileID = fopen(SummaryPath,'w');
assert(FileID >= 0,'Cannot create epsilon sweep summary.');
FileCleanup = onCleanup(@() fclose(FileID)); %#ok<NASGU>

fprintf(FileID,'# Control-task DAC epsilon sweep\n\n');
fprintf(FileID,['Fixed settings: tracking criterion, tracking tolerance ', ...
    '= 0.85, Kappa_P = 5, 400 inducing points, 30 s.\n\n']);
fprintf(FileID,['| Epsilon | Final control tracking error | Total broadcasts | ', ...
    'First 5 s trigger mean | Last 5 s trigger mean | ', ...
    'Mean DAC iterations | Zero-iteration updates |\n']);
fprintf(FileID,'|---:|---:|---:|---:|---:|---:|---:|\n');
for EpsilonNr = 1:numel(EpsilonSet)
    fprintf(FileID,'| %g | %.6g | %.0f | %.6g | %.6g | %.6g | %.0f |\n', ...
        EpsilonSet(EpsilonNr),FinalControlTrackingError(EpsilonNr), ...
        TotalBroadcasts(EpsilonNr),EarlyTriggerMean(EpsilonNr), ...
        LateTriggerMean(EpsilonNr),MeanDACIterations(EpsilonNr), ...
        ZeroIterationUpdates(EpsilonNr));
end

fprintf('\nDAC epsilon sweep\n');
fprintf(['Epsilon  FinalControlError  TotalBroadcasts  EarlyTrigger  ', ...
    'LateTrigger  MeanIterations  ZeroIterations\n']);
for EpsilonNr = 1:numel(EpsilonSet)
    fprintf('%7g  %17.6g  %15.0f  %12.6g  %11.6g  %14.6g  %14.0f\n', ...
        EpsilonSet(EpsilonNr),FinalControlTrackingError(EpsilonNr), ...
        TotalBroadcasts(EpsilonNr),EarlyTriggerMean(EpsilonNr), ...
        LateTriggerMean(EpsilonNr),MeanDACIterations(EpsilonNr), ...
        ZeroIterationUpdates(EpsilonNr));
end
fprintf('Markdown: %s\n',SummaryPath);

function restore_environment(Names,Values)
for EnvironmentNr = 1:numel(Names)
    setenv(Names{EnvironmentNr},Values{EnvironmentNr});
end
end
