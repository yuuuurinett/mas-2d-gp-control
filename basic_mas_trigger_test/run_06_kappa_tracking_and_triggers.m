%% Basic MAS: per-point frozen Zeta diffusion diagnostic
% Same as run_06_kappa_tracking_and_triggers.m, EXCEPT the inner loop no
% longer uses a global tolerance to decide when to stop, and no longer
% pushes Zeta uniformly across every inducing point.
%
% Equation (17) -- the actual broadcast trigger -- is completely
% unchanged. The only change is WHERE the Zeta diffusion step is allowed
% to act: at each fixed iteration, an (agent, point) pair is only pushed
% if its neighbors' last-broadcast values (XiHat) still disagree by more
% than ZetaWorkTolerance -- a near-exact-equality threshold, NOT the
% much looser DACTriggerOffset that equation (17) uses to decide
% broadcasts. Pairs that have already settled to (near) exact agreement
% stop being nudged, instead of every point being force-pushed for a
% fixed number of steps regardless of whether it has already converged.
%
% Diagnostic question: does the number of "still active" (agent,point)
% pairs decay over online-update-time on its own, now that convergence
% is judged per point instead of by one global before/after tolerance?

close all; clear; clc;

ScriptFolder = fileparts(mfilename('fullpath'));
ProjectRoot = fileparts(ScriptFolder);
addpath(ProjectRoot);
addpath(fullfile(ProjectRoot,'GP_inducingpoint'));

%% Settings
rng(0);
AgentQuantity = 6;
KappaP = 2;
TimeStep = 0.01;
SimulationEndTime = 40;
OnlineUpdateInterval = 0.1;
OnlineUpdateStep = round(OnlineUpdateInterval/TimeStep);
Time = 0:TimeStep:SimulationEndTime;
OnlineUpdateTime = Time(1:OnlineUpdateStep:end);

DiagnosticTolerance = 0.85;    % reporting only, does not stop the loop
DACTriggerOffset = 0.04;       % equation (17), unchanged
FixedDACSteps = 10;            % every update runs exactly this many steps
ZetaWorkTolerance = 3e-3;      % near-exact-equality freeze threshold
NumInducingPoints = 40;
MaxTrainingPoints = numel(OnlineUpdateTime)+10;

PlotEndTime = 30;

%% Undirected ring graph
Adjacency = zeros(AgentQuantity);
for AgentNr = 1:AgentQuantity
    NeighborNr = mod(AgentNr,AgentQuantity)+1;
    Adjacency(AgentNr,NeighborNr) = 1;
    Adjacency(NeighborNr,AgentNr) = 1;
end
Degree = sum(Adjacency,2);
Laplacian = diag(Degree)-Adjacency;

%% Simple MAS and LocalGP objects
AgentState = zeros(AgentQuantity,numel(Time));
AgentState(:,1) = rand(AgentQuantity,1);

StateDimension = 1;
OutputDimension = 2;
LocalGPSet = cell(AgentQuantity,1);
for AgentNr = 1:AgentQuantity
    LocalGPSet{AgentNr} = LocalGP_MultiOutput( ...
        StateDimension,OutputDimension,MaxTrainingPoints,0.01,1,0.5);
    LocalGPSet{AgentNr}.tau = 1e-8;
    LocalGPSet{AgentNr}.delta = 0.1;
end
InducingPointCoordinates = linspace(-2.5,3.5,NumInducingPoints);

%% Stored values
UpdateQuantity = numel(OnlineUpdateTime);
TrackingErrorBeforeDAC = nan(AgentQuantity,UpdateQuantity);
TrackingErrorAfterDAC = nan(AgentQuantity,UpdateQuantity);
AbsolutePChange = nan(1,UpdateQuantity);
TriggerCountPerAgent = zeros(AgentQuantity,UpdateQuantity);
% How many (agent,point) pairs were still "active" (being pushed by
% Zeta) at the LAST of the FixedDACSteps iterations in this update --
% this is the key new diagnostic.
ActivePairCountFinalStep = zeros(1,UpdateQuantity);
% Same, but averaged across all FixedDACSteps iterations in this update.
ActivePairCountMeanStep = zeros(1,UpdateQuantity);

CurrentP = [];
PreviousP = [];
Zeta = [];
XiHat = [];
Xi = [];
OnlineUpdateNr = 0;

%% Simulation
for TimeIndex = 1:numel(Time)
    IsOnlineUpdate = mod(TimeIndex-1,OnlineUpdateStep) == 0;

    if IsOnlineUpdate
        OnlineUpdateNr = OnlineUpdateNr+1;

        % 1. Every agent adds one new measurement.
        for AgentNr = 1:AgentQuantity
            CurrentState = AgentState(AgentNr,TimeIndex);
            TrueOutput = 0.2*[ ...
                sin(CurrentState)+0.15*CurrentState^2; ...
                cos(0.7*CurrentState)-1];
            MeasuredOutput = TrueOutput ...
                + 0.001*randn(OutputDimension,1);
            LocalGPSet{AgentNr}.addPoint(CurrentState,MeasuredOutput);
        end

        % 2. Recompute the reference signal P.
        [CurrentP,~] = gp_masked_aggregation_init( ...
            LocalGPSet,AgentQuantity,NumInducingPoints, ...
            InducingPointCoordinates,'poe');

        if ~isempty(PreviousP)
            PChange = CurrentP-PreviousP;
            PChangeNorm = sqrt(sum(PChange.^2,1));
            AbsolutePChange(OnlineUpdateNr) = max(PChangeNorm(:));
        end
        PreviousP = CurrentP;

        if isempty(Zeta)
            Zeta = zeros(size(CurrentP));
            XiHat = CurrentP;
        end
        Xi = CurrentP-Zeta;
        PAverage = mean(CurrentP,2);

        for AgentNr = 1:AgentQuantity
            Difference = Xi(:,AgentNr,:)-PAverage;
            PointNorm = sqrt(sum(Difference.^2,1));
            TrackingErrorBeforeDAC(AgentNr,OnlineUpdateNr) = ...
                max(PointNorm(:));
        end

        % 3. Fixed-step inner loop, per-(agent,point) Zeta freezing.
        % Equation (17) itself -- the broadcast decision below -- is
        % completely unchanged from run_06_kappa_tracking_and_triggers.m.
        ActivePairCountSet = zeros(1,FixedDACSteps);
        for IterationNr = 1:FixedDACSteps

            % Decide, for every (agent,point) pair, whether its
            % neighbors' last-broadcast values (XiHat) still disagree by
            % more than the near-exact-equality tolerance. Pairs that
            % have already settled are frozen for this step -- they get
            % no further Zeta push, so they cannot drift and re-trigger
            % on noise alone.
            ActiveMask = false(AgentQuantity,NumInducingPoints);
            for AgentNr = 1:AgentQuantity
                NeighborIndices = find(Adjacency(AgentNr,:) > 0);
                NeighborDisagreementSquared = zeros(1,NumInducingPoints);
                for NeighborNr = NeighborIndices
                    NeighborDifference = XiHat(:,AgentNr,:) ...
                        - XiHat(:,NeighborNr,:);
                    NeighborDisagreementSquared = ...
                        NeighborDisagreementSquared ...
                        + Adjacency(AgentNr,NeighborNr)*reshape( ...
                        sum(NeighborDifference.^2,1),1,NumInducingPoints);
                end
                ActiveMask(AgentNr,:) = ...
                    NeighborDisagreementSquared > ZetaWorkTolerance^2;
            end
            ActivePairCountSet(IterationNr) = nnz(ActiveMask);

            XiHatAgentFirst = permute(XiHat,[2,1,3]);
            XiHatFlat = reshape(XiHatAgentFirst,AgentQuantity,[]);
            LaplacianXiHatFlat = Laplacian*XiHatFlat;
            LaplacianXiHat = permute(reshape( ...
                LaplacianXiHatFlat,AgentQuantity, ...
                size(XiHat,1),NumInducingPoints),[2,1,3]);

            % Only push Zeta where ActiveMask is true; frozen pairs keep
            % their existing Zeta value unchanged this iteration.
            ActiveMask3D = reshape(ActiveMask,1,AgentQuantity,NumInducingPoints);
            Zeta = Zeta + TimeStep*KappaP*LaplacianXiHat.*ActiveMask3D;
            Xi = CurrentP-Zeta;

            for AgentNr = 1:AgentQuantity
                NeighborIndices = find(Adjacency(AgentNr,:) > 0);

                BroadcastError = XiHat(:,AgentNr,:)-Xi(:,AgentNr,:);
                BroadcastErrorSquared = reshape( ...
                    sum(BroadcastError.^2,1),1,NumInducingPoints);

                NeighborDisagreementSquared = zeros(1,NumInducingPoints);
                for NeighborNr = NeighborIndices
                    NeighborDifference = XiHat(:,AgentNr,:) ...
                        - XiHat(:,NeighborNr,:);
                    NeighborDisagreementSquared = ...
                        NeighborDisagreementSquared ...
                        + Adjacency(AgentNr,NeighborNr)*reshape( ...
                        sum(NeighborDifference.^2,1),1,NumInducingPoints);
                end

                CommunicationThreshold = ...
                    NeighborDisagreementSquared/(4*Degree(AgentNr)) ...
                    + DACTriggerOffset^2;
                TriggeredPoint = ...
                    BroadcastErrorSquared > CommunicationThreshold;

                if OnlineUpdateNr == 1 && IterationNr == 1
                    TriggeredPoint(:) = true;
                end

                if any(TriggeredPoint)
                    XiHat(:,AgentNr,TriggeredPoint) = ...
                        Xi(:,AgentNr,TriggeredPoint);
                    TriggerCountPerAgent(AgentNr,OnlineUpdateNr) = ...
                        TriggerCountPerAgent(AgentNr,OnlineUpdateNr) ...
                        + sum(TriggeredPoint);
                end
            end
        end

        ActivePairCountFinalStep(OnlineUpdateNr) = ActivePairCountSet(end);
        ActivePairCountMeanStep(OnlineUpdateNr) = mean(ActivePairCountSet);

        for AgentNr = 1:AgentQuantity
            Difference = Xi(:,AgentNr,:)-PAverage;
            PointNorm = sqrt(sum(Difference.^2,1));
            TrackingErrorAfterDAC(AgentNr,OnlineUpdateNr) = ...
                max(PointNorm(:));
        end
    end

    if TimeIndex < numel(Time)
        ControlInput = -KappaP*Laplacian*AgentState(:,TimeIndex);
        AgentState(:,TimeIndex+1) = AgentState(:,TimeIndex) ...
            + TimeStep*ControlInput;
    end
end

%% Result structure
Result = struct();
Result.KappaP = KappaP;
Result.Time = Time;
Result.OnlineUpdateTime = OnlineUpdateTime;
Result.FixedDACSteps = FixedDACSteps;
Result.ZetaWorkTolerance = ZetaWorkTolerance;
Result.DACTriggerOffset = DACTriggerOffset;
Result.NumInducingPoints = NumInducingPoints;
Result.TrackingErrorBeforeDAC = TrackingErrorBeforeDAC;
Result.TrackingErrorAfterDAC = TrackingErrorAfterDAC;
Result.AbsolutePChange = AbsolutePChange;
Result.TriggerCountPerAgent = TriggerCountPerAgent;
Result.TriggerCountPerAgentPerPoint = TriggerCountPerAgent/NumInducingPoints;
Result.ActivePairCountFinalStep = ActivePairCountFinalStep;
Result.ActivePairCountMeanStep = ActivePairCountMeanStep;

AgentColors = lines(AgentQuantity);
ResultTag = sprintf('kappa_%g_M%d_T%g_perpointfreeze', ...
    KappaP,NumInducingPoints,SimulationEndTime);

%% Figure 1: the key new diagnostic -- active (agent,point) pair count
ActiveFigure = figure('Color','w','Position',[80,60,1150,520]);
ActiveAxes = axes(ActiveFigure);
hold(ActiveAxes,'on');
plot(ActiveAxes,OnlineUpdateTime,ActivePairCountFinalStep,'-', ...
    'LineWidth',1.0,'DisplayName','Active pairs, last of 10 steps');
plot(ActiveAxes,OnlineUpdateTime,ActivePairCountMeanStep,'-', ...
    'LineWidth',1.0,'DisplayName','Active pairs, averaged over 10 steps');
MaxPossiblePairs = AgentQuantity*NumInducingPoints;
yline(ActiveAxes,MaxPossiblePairs,'--k', ...
    sprintf('All %d pairs active',MaxPossiblePairs),'HandleVisibility','off');
grid(ActiveAxes,'on'); box(ActiveAxes,'on');
xlim(ActiveAxes,[0,PlotEndTime]);
xlabel(ActiveAxes,'Time (s)');
ylabel(ActiveAxes,'Active (agent,point) pairs');
title(ActiveAxes, ...
    'Per-point freeze diagnostic: how many pairs still need pushing');
legend(ActiveAxes,'Location','northeast');

ActiveFigurePath = fullfile(ScriptFolder, ...
    ['result_06_active_pairs_',ResultTag,'_view0to',num2str(PlotEndTime),'s.png']);
exportgraphics(ActiveFigure,ActiveFigurePath,'Resolution',180);

%% Figure 2: P change
PChangeFigure = figure('Color','w','Position',[100,90,1150,500]);
PChangeAxes = axes(PChangeFigure);
ValidPChange = ~isnan(AbsolutePChange);
semilogy(PChangeAxes,OnlineUpdateTime(ValidPChange), ...
    AbsolutePChange(ValidPChange),'-','LineWidth',1.0);
grid(PChangeAxes,'on'); box(PChangeAxes,'on');
xlim(PChangeAxes,[0,PlotEndTime]);
xlabel(PChangeAxes,'Time (s)'); ylabel(PChangeAxes,'Absolute P change');
title(PChangeAxes,'Maximum absolute reference-signal change at each update');

PChangeFigurePath = fullfile(ScriptFolder, ...
    ['result_06_p_change_',ResultTag,'_view0to',num2str(PlotEndTime),'s.png']);
exportgraphics(PChangeFigure,PChangeFigurePath,'Resolution',180);

%% Figure 3: broadcasts in every 0.1 s update
CountFigure = figure('Color','w','Position',[100,80,1250,540]);
CountAxes = axes(CountFigure);
hold(CountAxes,'on');
for AgentNr = 1:AgentQuantity
    plot(CountAxes,OnlineUpdateTime, ...
        TriggerCountPerAgent(AgentNr,:)/NumInducingPoints,'-', ...
        'Color',AgentColors(AgentNr,:),'LineWidth',0.8, ...
        'DisplayName',sprintf('Agent %d',AgentNr));
end
grid(CountAxes,'on'); box(CountAxes,'on');
xlim(CountAxes,[0,PlotEndTime]);
xlabel(CountAxes,'Time (s)');
ylabel(CountAxes,'Broadcasts per inducing point');
title(CountAxes,'DAC broadcasts at every 0.1 s update (per-point freeze)');
legend(CountAxes,'Location','eastoutside');

CountFigurePath = fullfile(ScriptFolder, ...
    ['result_06_trigger_count_0p1s_',ResultTag,'_view0to',num2str(PlotEndTime),'s.png']);
exportgraphics(CountFigure,CountFigurePath,'Resolution',180);

%% Save data and Markdown summary
DataPath = fullfile(ScriptFolder,['result_06_',ResultTag,'.mat']);
save(DataPath,'Result');

MarkdownPath = fullfile(ScriptFolder,['result_06_',ResultTag,'_summary.md']);
FileID = fopen(MarkdownPath,'w');
assert(FileID >= 0,'Cannot create Markdown summary.');
CleanupFile = onCleanup(@() fclose(FileID));

fprintf(FileID,'# Per-point frozen Zeta diffusion diagnostic\n\n');
fprintf(FileID,'- Kappa_P: %g\n',KappaP);
fprintf(FileID,'- Fixed DAC steps per update: %d (no tolerance-based exit)\n', ...
    FixedDACSteps);
fprintf(FileID,'- DAC trigger offset (equation 17, unchanged): %g\n', ...
    DACTriggerOffset);
fprintf(FileID,'- Zeta work tolerance (freeze threshold): %g\n', ...
    ZetaWorkTolerance);
fprintf(FileID,'- Inducing points: %d\n', NumInducingPoints);
fprintf(FileID,'- Max possible active pairs: %d (AgentQuantity x NumInducingPoints)\n\n', ...
    AgentQuantity*NumInducingPoints);

fprintf(FileID,['| Interval (s) | Mean |ΔP| | Broadcasts per agent per point | ', ...
    'Mean active pairs (final step) | Mean active pairs (avg over steps) |\n']);
fprintf(FileID,'|---:|---:|---:|---:|---:|\n');
for WindowStart = 0:10:(SimulationEndTime-10)
    WindowMask = OnlineUpdateTime >= WindowStart ...
        & OnlineUpdateTime < WindowStart+10;
    MeanPChange = mean(AbsolutePChange(WindowMask),'omitnan');
    MeanBroadcast = sum(TriggerCountPerAgent(:,WindowMask),'all') ...
        /(AgentQuantity*NumInducingPoints*sum(WindowMask));
    MeanActiveFinal = mean(ActivePairCountFinalStep(WindowMask));
    MeanActiveAvg = mean(ActivePairCountMeanStep(WindowMask));
    fprintf(FileID,'| %g-%g | %.6g | %.6g | %.4g | %.4g |\n', ...
        WindowStart,WindowStart+10,MeanPChange,MeanBroadcast, ...
        MeanActiveFinal,MeanActiveAvg);
end

clear CleanupFile;

fprintf('Finished per-point freeze diagnostic.\n');
fprintf('Active-pairs figure: %s\n',ActiveFigurePath);
fprintf('P-change figure:     %s\n',PChangeFigurePath);
fprintf('Trigger count figure:%s\n',CountFigurePath);
fprintf('Markdown:            %s\n',MarkdownPath);
fprintf('Data:                %s\n',DataPath);