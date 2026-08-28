%% Step 2: tune only Kappa_P
% Fixed: tracking criterion, tracking tolerance = 0.85, epsilon = 0.005.

clear; clc;

ProjectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(ProjectRoot));

ResultFolder = fullfile(ProjectRoot,'control_test','result','Diagnostics');
if ~isfolder(ResultFolder)
    mkdir(ResultFolder);
end

KappaSet = [1,3,4,5];
ResultFileSet = { ...
    'control_task_dac_graft_tracking_tau085_30.mat'; ...
    'control_task_kappa3_tau085_eps0005_30.mat'; ...
    'control_task_kappa4_tau085_eps0005_30.mat'; ...
    'control_task_dac_graft_tracking_tau085_kappa5_30.mat'};

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
setenv('CONTROL_DAC_TRIGGER_EPSILON','0.005');
setenv('CONTROL_DAC_TRACKING_TOL','0.85');
setenv('CONTROL_DAC_CONSENSUS_TOL','0.005');
setenv('CONTROL_DAC_TERMINATION_CRITERION','tracking');

for KappaNr = 1:numel(KappaSet)
    ResultPath = fullfile(ResultFolder,ResultFileSet{KappaNr});
    if isfile(ResultPath)
        fprintf('Reuse Kappa_P=%g: %s\n',KappaSet(KappaNr),ResultPath);
        continue;
    end

    setenv('CONTROL_DAC_KAPPA_P',num2str(KappaSet(KappaNr),16));
    [~,ResultBaseName] = fileparts(ResultFileSet{KappaNr});
    fprintf('Run Kappa_P=%g ...\n',KappaSet(KappaNr));
    run_simulation_inducing_point( ...
        'poe',ResultFolder,ResultBaseName,true,0.2,0.1,0,42);
end

%% Requested metrics plus mean DAC iterations
FinalControlTrackingError = zeros(size(KappaSet));
TotalBroadcasts = zeros(size(KappaSet));
EarlyTriggerMean = zeros(size(KappaSet));
LateTriggerMean = zeros(size(KappaSet));
ZeroIterationUpdates = zeros(size(KappaSet));
MeanDACIterations = zeros(size(KappaSet));

for KappaNr = 1:numel(KappaSet)
    ResultPath = fullfile(ResultFolder,ResultFileSet{KappaNr});
    Data = load(ResultPath,'TrackingError_vector', ...
        'dac_event_count_per_agent','DACBroadcastCountPerSecond', ...
        'dac_iteration_count_set','consensus_communication_update_set');

    FinalControlTrackingError(KappaNr) = Data.TrackingError_vector(end);
    TotalBroadcasts(KappaNr) = sum(Data.dac_event_count_per_agent);
    EarlyTriggerMean(KappaNr) = ...
        mean(Data.DACBroadcastCountPerSecond(1:5));
    LateTriggerMean(KappaNr) = ...
        mean(Data.DACBroadcastCountPerSecond(end-4:end));

    OnlineUpdateMask = ...
        Data.consensus_communication_update_set(1:end-1) > 0;
    UpdateIterationCount = ...
        Data.dac_iteration_count_set(OnlineUpdateMask);
    ZeroIterationUpdates(KappaNr) = sum(UpdateIterationCount == 0);
    MeanDACIterations(KappaNr) = mean(UpdateIterationCount);
end

SummaryPath = fullfile(ResultFolder,'control_kappa_sweep_summary.md');
FileID = fopen(SummaryPath,'w');
assert(FileID >= 0,'Cannot create Kappa_P sweep summary.');
FileCleanup = onCleanup(@() fclose(FileID)); %#ok<NASGU>

fprintf(FileID,'# Control-task Kappa_P sweep\n\n');
fprintf(FileID,['Fixed settings: tracking criterion, tracking tolerance ', ...
    '= 0.85, epsilon = 0.005, 400 inducing points, 30 s.\n\n']);
fprintf(FileID,['| Kappa_P | Final control tracking error | Total broadcasts | ', ...
    'First 5 s trigger mean | Last 5 s trigger mean | ', ...
    'Mean DAC iterations | Zero-iteration updates |\n']);
fprintf(FileID,'|---:|---:|---:|---:|---:|---:|---:|\n');
for KappaNr = 1:numel(KappaSet)
    fprintf(FileID,'| %g | %.6g | %.0f | %.6g | %.6g | %.6g | %.0f |\n', ...
        KappaSet(KappaNr),FinalControlTrackingError(KappaNr), ...
        TotalBroadcasts(KappaNr),EarlyTriggerMean(KappaNr), ...
        LateTriggerMean(KappaNr),MeanDACIterations(KappaNr), ...
        ZeroIterationUpdates(KappaNr));
end

fprintf('\nKappa_P sweep\n');
fprintf(['Kappa  FinalControlError  TotalBroadcasts  EarlyTrigger  ', ...
    'LateTrigger  MeanIterations  ZeroIterations\n']);
for KappaNr = 1:numel(KappaSet)
    fprintf('%5g  %17.6g  %15.0f  %12.6g  %11.6g  %14.6g  %14.0f\n', ...
        KappaSet(KappaNr),FinalControlTrackingError(KappaNr), ...
        TotalBroadcasts(KappaNr),EarlyTriggerMean(KappaNr), ...
        LateTriggerMean(KappaNr),MeanDACIterations(KappaNr), ...
        ZeroIterationUpdates(KappaNr));
end
fprintf('Markdown: %s\n',SummaryPath);

function restore_environment(Names,Values)
for EnvironmentNr = 1:numel(Names)
    setenv(Names{EnvironmentNr},Values{EnvironmentNr});
end
end
