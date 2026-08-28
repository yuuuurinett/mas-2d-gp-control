%% Plot the control-task result using the same quantities as simple MAS

clear; clc;

ScriptFolder = fileparts(mfilename('fullpath'));
ResultFolder = fullfile(ScriptFolder,'result','Diagnostics');
ResultPath = fullfile(ResultFolder, ...
    'control_task_eps005_tau085_kappa5_30.mat');

Result = load(ResultPath,'t_set','AbsolutePChangeRMS', ...
    'dac_broadcast_count_set','NumInducingPoints','DACTriggerEpsilon', ...
    'DACTrackingTolerance','Kappa_P');

AgentQuantity = size(Result.dac_broadcast_count_set,1);
AgentColors = lines(AgentQuantity);
PhysicalTime = Result.t_set(1:end-1);

% P change exists only at the 0.1 s online updates.
PUpdateMask = isfinite(Result.AbsolutePChangeRMS(1:end-1));
PUpdateTime = PhysicalTime(PUpdateMask);
PChange = Result.AbsolutePChangeRMS(PUpdateMask);

% Each saved entry is already the average number of broadcasts over all
% inducing points for one agent during one 0.1 s online update.
BroadcastPerAgentPointAtUpdate = Result.dac_broadcast_count_set;
OnlineUpdateMask = any(BroadcastPerAgentPointAtUpdate > 0,1);
OnlineUpdateTime = PhysicalTime(OnlineUpdateMask);
BroadcastAtUpdate = BroadcastPerAgentPointAtUpdate(:,OnlineUpdateMask);

SimulationEndTime = Result.t_set(end);
SecondEdges = 0:1:SimulationEndTime;
SecondCenters = SecondEdges(1:end-1)+0.5;
BroadcastPerSecond = zeros(AgentQuantity,numel(SecondCenters));
for SecondNr = 1:numel(SecondCenters)
    WindowMask = PhysicalTime >= SecondEdges(SecondNr) & ...
        PhysicalTime < SecondEdges(SecondNr+1);
    BroadcastPerSecond(:,SecondNr) = ...
        sum(BroadcastPerAgentPointAtUpdate(:,WindowMask),2);
end

%% Figure 1: absolute P change at every online update
PChangeFigure = figure('Color','w','Position',[100,100,1150,500]);
semilogy(PUpdateTime,PChange,'-','LineWidth',1.1);
grid on; box on; xlim([0,SimulationEndTime]);
xlabel('Time (s)'); ylabel('RMS of P(k)-P(k-1)');
title('Control task: absolute reference-signal change at each 0.1 s update');
PChangePath = fullfile(ResultFolder,'control_like_simple_p_change.png');
exportgraphics(PChangeFigure,PChangePath,'Resolution',180);
close(PChangeFigure);

%% Figure 2: broadcasts at every 0.1 s update, six agents
UpdateFigure = figure('Color','w','Position',[100,100,1250,540]);
hold on;
for AgentNr = 1:AgentQuantity
    plot(OnlineUpdateTime,BroadcastAtUpdate(AgentNr,:),'-', ...
        'Color',AgentColors(AgentNr,:),'LineWidth',0.9, ...
        'DisplayName',sprintf('Agent %d',AgentNr));
end
grid on; box on; xlim([0,SimulationEndTime]);
xlabel('Time (s)');
ylabel('Broadcasts per inducing point');
title('Control task: DAC broadcasts at every 0.1 s online update');
legend('Location','eastoutside');
UpdatePath = fullfile(ResultFolder,'control_like_simple_trigger_count_0p1s.png');
exportgraphics(UpdateFigure,UpdatePath,'Resolution',180);
close(UpdateFigure);

%% Figure 3: simple-MAS-style star plot of one-second averages
StarFigure = figure('Color','w','Position',[80,80,1250,520]);
hold on;
for AgentNr = 1:AgentQuantity
    StarQuantity = round(BroadcastPerSecond(AgentNr,:));
    for SecondNr = 1:numel(SecondCenters)
        if StarQuantity(SecondNr) > 0
            StarTime = linspace(SecondEdges(SecondNr)+0.08, ...
                SecondEdges(SecondNr+1)-0.08,StarQuantity(SecondNr));
            plot(StarTime,AgentNr*ones(size(StarTime)),'*', ...
                'LineStyle','none','Color',AgentColors(AgentNr,:), ...
                'MarkerSize',4);
        end
    end
end
grid on; box on;
xlim([0,SimulationEndTime]); ylim([0.5,AgentQuantity+0.5]);
yticks(1:AgentQuantity);
xlabel('Time (s)'); ylabel('Agent');
title(['Average DAC broadcasts per agent per inducing point ', ...
    '(one star = one average broadcast)']);
StarPath = fullfile(ResultFolder,'control_like_simple_trigger_stars.png');
exportgraphics(StarFigure,StarPath,'Resolution',180);
close(StarFigure);

%% Word-friendly numerical summary
EarlyTimeMask = PUpdateTime < 5;
LateTimeMask = PUpdateTime >= SimulationEndTime-5;
EarlyPChange = mean(PChange(EarlyTimeMask));
LatePChange = mean(PChange(LateTimeMask));
EarlyTrigger = mean(BroadcastPerSecond(:,1:5),'all');
LateTrigger = mean(BroadcastPerSecond(:,end-4:end),'all');

SummaryPath = fullfile(ResultFolder, ...
    'control_like_simple_p_and_trigger_summary.md');
FileID = fopen(SummaryPath,'w');
assert(FileID >= 0,'Cannot create control-task plot summary.');
FileCleanup = onCleanup(@() fclose(FileID)); %#ok<NASGU>

fprintf(FileID,'# Control task: P change and DAC triggers\n\n');
fprintf(FileID,['Settings: tracking tolerance = %.3g, Kappa_P = %.3g, ', ...
    'epsilon = %.3g, inducing points = %d.\n\n'], ...
    Result.DACTrackingTolerance,Result.Kappa_P,Result.DACTriggerEpsilon, ...
    Result.NumInducingPoints);
fprintf(FileID,'| Quantity | First 5 s mean | Last 5 s mean | Change |\n');
fprintf(FileID,'|---|---:|---:|---:|\n');
fprintf(FileID,'| Absolute P-change RMS | %.6g | %.6g | %.2f%% |\n', ...
    EarlyPChange,LatePChange,100*(LatePChange/EarlyPChange-1));
fprintf(FileID,['| DAC broadcasts per agent per inducing point per second ', ...
    '| %.6g | %.6g | %.2f%% |\n'], ...
    EarlyTrigger,LateTrigger,100*(LateTrigger/EarlyTrigger-1));

fprintf('\nFirst 5 s -> last 5 s\n');
fprintf('Absolute P-change RMS: %.6g -> %.6g (%+.2f%%)\n', ...
    EarlyPChange,LatePChange,100*(LatePChange/EarlyPChange-1));
fprintf(['DAC broadcasts / agent / inducing point / second: ', ...
    '%.6g -> %.6g (%+.2f%%)\n'], ...
    EarlyTrigger,LateTrigger,100*(LateTrigger/EarlyTrigger-1));
fprintf('P-change plot: %s\n',PChangePath);
fprintf('0.1 s trigger plot: %s\n',UpdatePath);
fprintf('Star plot: %s\n',StarPath);
fprintf('Markdown: %s\n',SummaryPath);
