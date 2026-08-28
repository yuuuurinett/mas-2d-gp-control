function run_ip_reachable_uniform_case()
%RUN_IP_REACHABLE_UNIFORM_CASE Uniform M=400 points in reachable-state box.

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result','ip_reachable_uniform');
if ~isfolder(result_folder), mkdir(result_folder); end

baseline_file = fullfile(repo_root,'result','ip_agent_kappa_sweep', ...
    'poe_M400_R1_eps_0p1_kappa_1_T30.mat');
baseline = load(baseline_file,'x_all_set');
state_samples = reshape(baseline.x_all_set,4,[]);
state_min = min(state_samples,[],2);
state_max = max(state_samples,[],2);
state_center = (state_min+state_max)/2;
state_half_width = 1.10*(state_max-state_min)/2;
sampling_min = state_center-state_half_width;
sampling_max = state_center+state_half_width;

rng(42,'twister');
M = 400;
InducingPoints_Coordinates = sampling_min + ...
    (sampling_max-sampling_min).*rand(4,M);
point_file = fullfile(result_folder,'reachable_uniform_M400.mat');
save(point_file,'InducingPoints_Coordinates','state_min','state_max', ...
    'sampling_min','sampling_max');

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_DAC_STEPS_PER_PHYSICAL_STEP','CONTROL_IP_PROJECTION_DIAGNOSTICS', ...
    'CONTROL_IP_INDUCING_POINT_FILE','CONTROL_DAC_TRIGGER_EPSILON', ...
    'CONTROL_DAC_KAPPA_P'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment(env_names,old_values)); %#ok<NASGU>
setenv(env_names{1},'30');
setenv(env_names{2},'400');
setenv(env_names{3},'1');
setenv(env_names{4},'1');
setenv(env_names{5},point_file);
setenv(env_names{6},'0.1');
setenv(env_names{7},'1');

run_simulation_inducing_point('poe',result_folder, ...
    'poe_M400_reachable_uniform_R1_eps_0p1_kappa_1_T30', ...
    false,[],[],[],42);
plot_comparison(baseline_file,fullfile(result_folder, ...
    'poe_M400_reachable_uniform_R1_eps_0p1_kappa_1_T30.mat'), ...
    result_folder);
end

function plot_comparison(baseline_file,reachable_file,result_folder)
base = load(baseline_file);
reach = load(reachable_file);
labels = {'Full-domain uniform';'Reachable-box uniform'};
sources = {base;reach};
MeanDirect = nan(2,1); MeanIPCEN = nan(2,1); MeanIPDAC = nan(2,1);
MeanTracking = nan(2,1); BroadcastRate = nan(2,1);
for idx = 1:2
    data = sources{idx};
    valid = isfinite(data.direct_prediction_error_norm_vector) & ...
        isfinite(data.ip_cen_prediction_error_norm_vector) & ...
        isfinite(data.prediction_error_norm_vector);
    MeanDirect(idx) = mean(data.direct_prediction_error_norm_vector(valid));
    MeanIPCEN(idx) = mean(data.ip_cen_prediction_error_norm_vector(valid));
    MeanIPDAC(idx) = mean(data.prediction_error_norm_vector(valid));
    MeanTracking(idx) = trapz(data.t_set,data.TrackingError_vector)/data.t_set(end);
    BroadcastRate(idx) = sum(data.dac_inner_trigger_instance_set,'all') / ...
        (size(data.dac_inner_trigger_instance_set,1)*data.t_set(end));
end
Summary = table(labels,MeanDirect,MeanIPCEN,MeanIPDAC,MeanTracking, ...
    BroadcastRate,'VariableNames',{'Sampling','MeanDirectError', ...
    'MeanIPCENError','MeanIPDACError','MeanTrackingError', ...
    'BroadcastsPerAgentSecond'});
writetable(Summary,fullfile(result_folder,'comparison.csv'));
save(fullfile(result_folder,'comparison.mat'),'Summary');

fig = figure('Color','w','Position',[100 100 1180 500]);
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(layout); hold(ax,'on');
plot(ax,reach.t_set,reach.direct_prediction_error_norm_vector, ...
    'LineWidth',1.4,'DisplayName','Direct PoE');
plot(ax,reach.t_set,reach.ip_cen_prediction_error_norm_vector, ...
    'LineWidth',1.4,'DisplayName','IP-CEN');
plot(ax,reach.t_set,reach.prediction_error_norm_vector, ...
    'LineWidth',1.4,'DisplayName','IP-DAC');
xlabel(ax,'Time (s)'); ylabel(ax,'Prediction error norm');
title(ax,'Reachable-box uniform M=400'); grid(ax,'on'); box(ax,'on'); legend(ax);

ax = nexttile(layout);
bar(ax,[MeanDirect,MeanIPCEN,MeanIPDAC]);
set(ax,'XTickLabel',labels); ylabel(ax,'Mean prediction error');
legend(ax,{'Direct PoE','IP-CEN','IP-DAC'},'Location','best');
title(ax,'Uniform sampling-domain comparison'); grid(ax,'on'); box(ax,'on');
exportgraphics(fig,fullfile(result_folder, ...
    'reachable_uniform_comparison.png'),'Resolution',230);
savefig(fig,fullfile(result_folder,'reachable_uniform_comparison.fig'));
disp(Summary);
end

function restore_environment(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
