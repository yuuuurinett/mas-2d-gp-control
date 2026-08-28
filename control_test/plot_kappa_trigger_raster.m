%% DAC trigger raster for the Kappa_P time-scale comparison
clear; close all; clc;

script_folder = fileparts(mfilename('fullpath'));
result_folder = fullfile(script_folder,'result','Diagnostics');

KappaValues = [1,5,10];
ResultFiles = { ...
    'persistent_dac_poe_continuous30.mat', ...
    'persistent_dac_poe_continuous30_kappa_05.mat', ...
    'persistent_dac_poe_continuous30_kappa_10.mat'};

AgentColors = lines(6);
figure('Color','w','Position',[80,40,1250,900]);

for CaseNr = 1:numel(KappaValues)
    Result = load(fullfile(result_folder,ResultFiles{CaseNr}), ...
        't_set','dac_broadcast_count_set');
    EventTime = Result.t_set(1:size(Result.dac_broadcast_count_set,2));
    EventOccurred = Result.dac_broadcast_count_set > 0;

    subplot(numel(KappaValues),1,CaseNr); hold on;
    LegendText = strings(1,size(EventOccurred,1));

    for AgentNr = 1:size(EventOccurred,1)
        EventIndices = find(EventOccurred(AgentNr,:));
        plot(EventTime(EventIndices), ...
            AgentNr*ones(size(EventIndices)), ...
            '*','Color',AgentColors(AgentNr,:), ...
            'MarkerSize',4,'LineWidth',0.8);
        LegendText(AgentNr) = sprintf('Agent %d, n=%d', ...
            AgentNr,numel(EventIndices));
    end

    grid on;
    xlim([EventTime(1),EventTime(end)]);
    ylim([0.5,size(EventOccurred,1)+0.5]);
    yticks(1:size(EventOccurred,1));
    ylabel('Agent');
    title(sprintf('POE-DAC, \\kappa_P = %g',KappaValues(CaseNr)));
    legend(LegendText,'Location','eastoutside','FontSize',8);

    if CaseNr == numel(KappaValues)
        xlabel('Time (s)');
    end
end

sgtitle('DAC broadcast trigger times under continuous online learning', ...
    'FontSize',14,'FontWeight','bold');

OutputPath = fullfile(result_folder, ...
    'kappa_comparison_trigger_raster.png');
exportgraphics(gcf,OutputPath,'Resolution',180);
fprintf('Figure saved: %s\n',OutputPath);
