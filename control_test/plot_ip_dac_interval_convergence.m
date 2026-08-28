function Stats = plot_ip_dac_interval_convergence(result_file)
%PLOT_IP_DAC_INTERVAL_CONVERGENCE Diagnose DAC inside each 0.1-s P interval.

if nargin < 1 || isempty(result_file)
    repo_root = fileparts(fileparts(mfilename('fullpath')));
    result_file = fullfile(repo_root,'result','ip_prediction_reference_comparison', ...
        'poe_M600_R10.mat');
end
data = load(result_file);
[result_folder,result_stem] = fileparts(result_file);

t = data.t_set(1:end-1);
trigger_count = double(data.dac_trigger_count_per_agent_point_set);
[agent_count,point_count,step_count] = size(trigger_count);
round_count = data.DACStepsPerPhysicalStep;
point_trigger_fraction = squeeze(sum(trigger_count,[1 2]))' / ...
    (agent_count*point_count*round_count);
agent_active = squeeze(any(trigger_count > 0,2));
agent_activation_fraction = mean(agent_active,1);
dac_error = data.dac_tracking_error_set;

online_steps = find(data.projection_update_set(1:step_count)>0);
if numel(online_steps) < 2
    error('At least two projection updates are required.');
end
interval_steps = round(median(diff(online_steps)));
phase_index = mod((0:step_count-1),interval_steps)+1;
phase_time = (0:interval_steps-1)*(data.t_set(2)-data.t_set(1));

early_mask = t < min(2,t(end)/2);
late_mask = t >= max(t(end)-2,t(end)/2);
early_point_phase = phase_average(point_trigger_fraction,phase_index,interval_steps,early_mask);
late_point_phase = phase_average(point_trigger_fraction,phase_index,interval_steps,late_mask);
early_error_phase = phase_geomean(dac_error,phase_index,interval_steps,early_mask);
late_error_phase = phase_geomean(dac_error,phase_index,interval_steps,late_mask);

period_count = ceil(t(end));
second_centers = (0:period_count-1)+0.5;
point_rate_second = nan(1,period_count);
agent_active_second = nan(1,period_count);
for second_i = 1:period_count
    mask = t >= second_i-1 & t < second_i;
    point_rate_second(second_i) = mean(point_trigger_fraction(mask));
    agent_active_second(second_i) = mean(agent_activation_fraction(mask));
end

fig = figure('Color','w','Position',[100 70 1180 760]);
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

rate_ax = nexttile(layout); hold(rate_ax,'on');
plot(rate_ax,second_centers,point_rate_second,'LineWidth',1.9);
xlabel(rate_ax,'Time (s)'); ylabel(rate_ax,'Broadcasts / agent / DAC round');
title(rate_ax,'Agent-level packaged broadcast activity'); grid(rate_ax,'on'); box(rate_ax,'on');

agent_ax = nexttile(layout); hold(agent_ax,'on');
plot(agent_ax,second_centers,agent_active_second,'LineWidth',1.9);
xlabel(agent_ax,'Time (s)'); ylabel(agent_ax,'Agent activation ratio');
ylim(agent_ax,[0 1.05]); title(agent_ax,'Physical steps with at least one broadcast');
grid(agent_ax,'on'); box(agent_ax,'on');

phase_trigger_ax = nexttile(layout); hold(phase_trigger_ax,'on');
plot(phase_trigger_ax,phase_time,early_point_phase,'LineWidth',1.9, ...
    'DisplayName','Early 0--2 s');
plot(phase_trigger_ax,phase_time,late_point_phase,'LineWidth',1.9, ...
    'DisplayName','Late final 2 s');
xlabel(phase_trigger_ax,'Time since P refresh (s)');
ylabel(phase_trigger_ax,'Broadcasts / agent / DAC round');
title(phase_trigger_ax,'Within one 0.1-s online interval');
legend(phase_trigger_ax,'Location','best'); grid(phase_trigger_ax,'on'); box(phase_trigger_ax,'on');

phase_error_ax = nexttile(layout); hold(phase_error_ax,'on');
semilogy(phase_error_ax,phase_time,early_error_phase,'LineWidth',1.9, ...
    'DisplayName','Early 0--2 s');
semilogy(phase_error_ax,phase_time,late_error_phase,'LineWidth',1.9, ...
    'DisplayName','Late final 2 s');
xlabel(phase_error_ax,'Time since P refresh (s)');
ylabel(phase_error_ax,'DAC average-tracking error');
title(phase_error_ax,'Consensus tracking inside interval');
legend(phase_error_ax,'Location','best'); grid(phase_error_ax,'on'); box(phase_error_ax,'on');

title(layout,sprintf('IP-DAC interval diagnostic: M=%d, %d DAC rounds / physical step', ...
    point_count,round_count));
output_base = fullfile(result_folder,[result_stem,'_interval_diagnostic']);
exportgraphics(fig,[output_base,'.png'],'Resolution',220);
savefig(fig,[output_base,'.fig']);

Window = {'Early 0-2 s';'Late final 2 s'};
MeanPointTriggerFractionPerRound = [mean(point_trigger_fraction(early_mask)); ...
    mean(point_trigger_fraction(late_mask))];
MeanAgentActivationRatio = [mean(agent_activation_fraction(early_mask)); ...
    mean(agent_activation_fraction(late_mask))];
MeanDACTrackingError = [mean(dac_error(early_mask),'omitnan'); ...
    mean(dac_error(late_mask),'omitnan')];
EndOfIntervalDACTrackingError = [early_error_phase(end);late_error_phase(end)];
Stats = table(Window,MeanPointTriggerFractionPerRound,MeanAgentActivationRatio, ...
    MeanDACTrackingError,EndOfIntervalDACTrackingError);
writetable(Stats,[output_base,'.csv']);
disp(Stats);
end

function output = phase_average(values,phase_index,phase_count,time_mask)
output = nan(1,phase_count);
for phase_i = 1:phase_count
    mask = phase_index == phase_i & time_mask;
    output(phase_i) = mean(values(mask),'omitnan');
end
end

function output = phase_geomean(values,phase_index,phase_count,time_mask)
output = nan(1,phase_count);
for phase_i = 1:phase_count
    mask = phase_index == phase_i & time_mask & isfinite(values) & values > 0;
    output(phase_i) = exp(mean(log(values(mask)),'omitnan'));
end
end
