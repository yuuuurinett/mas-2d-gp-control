%% Absolute P-update size versus DAC communication count
clear; close all; clc;

script_folder = fileparts(mfilename('fullpath'));
result_folder = fullfile(script_folder,'result','Diagnostics');
result_file = fullfile(result_folder, ...
    'persistent_dac_poe_continuous30.mat');
Result = load(result_file, ...
    't_set','projection_step_change_per_agent_set', ...
    'DACBroadcastCountPerSecond');

MeanAbsolutePChange = mean( ...
    Result.projection_step_change_per_agent_set,1,'omitnan');
MaxAbsolutePChange = max( ...
    Result.projection_step_change_per_agent_set,[],1);
ValidUpdate = isfinite(MeanAbsolutePChange);
UpdateTimes = Result.t_set(ValidUpdate);

figure('Color','w','Position',[100,80,1100,700]);

subplot(2,1,1);
semilogy(UpdateTimes,MeanAbsolutePChange(ValidUpdate), ...
    'o-','LineWidth',1.1,'MarkerSize',3, ...
    'DisplayName','Mean across agents');
hold on;
semilogy(UpdateTimes,MaxAbsolutePChange(ValidUpdate), ...
    '-','LineWidth',1.2, ...
    'DisplayName','Maximum agent');
grid on; xlim([0,Result.t_set(end)]);
xlabel('Time (s)'); ylabel('Absolute P-change RMS');
title('RMS of P(k)-P(k-1) at each 0.1 s online update');
legend('Location','northeast');

subplot(2,1,2);
WindowCentres = 0.5+(0:numel(Result.DACBroadcastCountPerSecond)-1);
bar(WindowCentres,Result.DACBroadcastCountPerSecond,1);
grid on; xlim([0,Result.t_set(end)]);
xlabel('One-second window centre');
ylabel('Total DAC broadcasts');
title('All-agent DAC broadcasts in each 1 s window');

sgtitle('Control task: absolute P update versus DAC communication', ...
    'FontSize',14,'FontWeight','bold');

OutputPath = fullfile(result_folder, ...
    'absolute_p_change_vs_dac_triggers.png');
exportgraphics(gcf,OutputPath,'Resolution',180);
fprintf('Figure saved: %s\n',OutputPath);
