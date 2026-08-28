function Summary = run_ip_inducing_location_ablation(do_run,m_value,dac_rounds)
%RUN_IP_INDUCING_LOCATION_ABLATION Uniform-domain versus trajectory IPs.
%
% The trajectory-box set is a diagnostic oracle constructed from the state
% range of the previously saved uniform-IP trajectory. It tests whether IP
% placement is the bottleneck; it is not a final online placement policy.

if nargin < 1 || isempty(do_run), do_run = true; end
if nargin < 2 || isempty(m_value), m_value = 600; end
if nargin < 3 || isempty(dac_rounds), dac_rounds = 10; end

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
baseline_folder = fullfile(repo_root,'result','ip_prediction_reference_comparison');
baseline_file = fullfile(baseline_folder,sprintf('poe_M%d_R%d.mat',m_value,dac_rounds));
result_folder = fullfile(repo_root,'result','ip_inducing_location_ablation');
if ~isfolder(result_folder), mkdir(result_folder); end

baseline = load(baseline_file);
trajectory_ip_file = fullfile(result_folder, ...
    sprintf('trajectory_box_inducing_points_M%d.mat',m_value));
InducingPoints_Coordinates = select_trajectory_box_points( ...
    baseline.x_all_set,m_value); %#ok<NASGU>
save(trajectory_ip_file,'InducingPoints_Coordinates');

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_DAC_STEPS_PER_PHYSICAL_STEP','CONTROL_IP_PROJECTION_DIAGNOSTICS', ...
    'CONTROL_IP_INDUCING_POINT_FILE'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment(env_names,old_values)); %#ok<NASGU>
setenv(env_names{1},num2str(baseline.t_set(end),'%.15g'));
setenv(env_names{2},num2str(m_value));
setenv(env_names{3},num2str(dac_rounds));
setenv(env_names{4},'1');
setenv(env_names{5},trajectory_ip_file);

trajectory_file = fullfile(result_folder, ...
    sprintf('trajectory_box_poe_M%d_R%d.mat',m_value,dac_rounds));
if do_run
    run_simulation_inducing_point('poe',result_folder, ...
        sprintf('trajectory_box_poe_M%d_R%d',m_value,dac_rounds), ...
        false,[],[],[],42);
end
trajectory = load(trajectory_file);

labels = {'Uniform-domain IP';'Trajectory-box IP'};
sets = {baseline;trajectory};
MeanDirectPoEError = nan(2,1);
MeanIPCENError = nan(2,1);
MeanIPDACError = nan(2,1);
FinalTrackingError = nan(2,1);
MeanPointTriggersPerAgentPointSecond = nan(2,1);
AgentActivationRatio = nan(2,1);
for idx = 1:2
    data = sets{idx};
    valid = isfinite(data.prediction_error_norm_vector);
    MeanDirectPoEError(idx) = mean(data.direct_prediction_error_norm_vector(valid));
    MeanIPCENError(idx) = mean(data.ip_cen_prediction_error_norm_vector(valid));
    MeanIPDACError(idx) = mean(data.prediction_error_norm_vector(valid));
    FinalTrackingError(idx) = data.TrackingError_vector(end);
    trigger_data = double(data.dac_trigger_count_per_agent_point_set);
    duration = data.t_set(end)-data.t_set(1);
    MeanPointTriggersPerAgentPointSecond(idx) = ...
        sum(trigger_data,'all')/(size(trigger_data,1)*size(trigger_data,2)*duration);
    agent_active = squeeze(any(trigger_data > 0,2));
    AgentActivationRatio(idx) = mean(agent_active,'all');
end

InducingPointPlacement = labels;
Summary = table(InducingPointPlacement,MeanDirectPoEError,MeanIPCENError, ...
    MeanIPDACError,FinalTrackingError,MeanPointTriggersPerAgentPointSecond, ...
    AgentActivationRatio);
writetable(Summary,fullfile(result_folder,sprintf('summary_M%d_R%d.csv',m_value,dac_rounds)));
save(fullfile(result_folder,sprintf('summary_M%d_R%d.mat',m_value,dac_rounds)), ...
    'Summary','m_value','dac_rounds');

plot_prediction_errors(baseline,trajectory,result_folder,m_value,dac_rounds);
plot_trigger_comparison(baseline,trajectory,result_folder,m_value,dac_rounds);
plot_trajectory_dynamics(trajectory,result_folder,m_value,dac_rounds);
disp(Summary);
end

function Z = select_trajectory_box_points(x_all_set,m_value)
x_dim = 4;
agent_count = size(x_all_set,1)/x_dim;
state_cube = reshape(x_all_set,x_dim,agent_count,[]);
candidates = reshape(state_cube,x_dim,[]);
candidates = candidates(:,all(isfinite(candidates),1));

% A trajectory itself is nearly a low-dimensional manifold; placing all 600
% points directly on it makes the squared-exponential kernel nearly singular.
% Instead, cover a padded full-dimensional box around the observed trajectory
% using deterministic stratified sampling in every coordinate.
lower = min(candidates,[],2);
upper = max(candidates,[],2);
center = 0.5*(lower+upper);
half_width = max(0.6*(upper-lower),0.15*ones(x_dim,1));
lower = center-half_width;
upper = center+half_width;

rng(20260813);
unit_points = zeros(x_dim,m_value);
for dimension_i = 1:x_dim
    strata = randperm(m_value);
    unit_points(dimension_i,:) = (strata-rand(1,m_value))/m_value;
end
Z = lower+(upper-lower).*unit_points;
end

function plot_prediction_errors(uniform_data,trajectory_data,result_folder,m_value,dac_rounds)
fig = figure('Color','w','Position',[100 80 1200 670]);
layout = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
plot_one_prediction_panel(nexttile(layout),uniform_data,'Uniform-domain inducing points');
plot_one_prediction_panel(nexttile(layout),trajectory_data,'Trajectory-box inducing points');
xlabel(layout,'Time (s)');
title(layout,sprintf('Effect of inducing-point placement (M=%d, DAC rounds=%d)', ...
    m_value,dac_rounds));
exportgraphics(fig,fullfile(result_folder, ...
    sprintf('prediction_location_comparison_M%d_R%d.png',m_value,dac_rounds)), ...
    'Resolution',220);
savefig(fig,fullfile(result_folder, ...
    sprintf('prediction_location_comparison_M%d_R%d.fig',m_value,dac_rounds)));
end

function plot_one_prediction_panel(ax,data,panel_title)
valid = isfinite(data.prediction_error_norm_vector);
t = data.t_set(valid);
hold(ax,'on');
plot(ax,t,data.direct_prediction_error_norm_vector(valid),'LineWidth',1.7, ...
    'DisplayName','Direct PoE');
plot(ax,t,data.ip_cen_prediction_error_norm_vector(valid),'LineWidth',1.7, ...
    'DisplayName','IP-CEN');
plot(ax,t,data.prediction_error_norm_vector(valid),'LineWidth',1.7, ...
    'DisplayName','IP-DAC');
ylabel(ax,'Prediction error norm'); title(ax,panel_title);
xlim(ax,[t(1),t(end)]); grid(ax,'on'); box(ax,'on');
legend(ax,'Location','best');
end

function plot_trigger_comparison(uniform_data,trajectory_data,result_folder,m_value,dac_rounds)
fig = figure('Color','w','Position',[120 100 1100 620]);
layout = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
rate_ax = nexttile(layout); hold(rate_ax,'on');
activation_ax = nexttile(layout); hold(activation_ax,'on');
plot_one_trigger_set(rate_ax,activation_ax,uniform_data,'Uniform-domain IP');
plot_one_trigger_set(rate_ax,activation_ax,trajectory_data,'Trajectory-box IP');
ylabel(rate_ax,'Point triggers / agent / point / s');
ylabel(activation_ax,'Agent activation ratio'); ylim(activation_ax,[0 1.05]);
xlabel(layout,'Time window (s)');
title(rate_ax,'Point-level event-trigger activity');
title(activation_ax,'Agent-level physical-step communication activity');
legend(rate_ax,'Location','best'); legend(activation_ax,'Location','best');
grid(rate_ax,'on'); box(rate_ax,'on'); grid(activation_ax,'on'); box(activation_ax,'on');
title(layout,sprintf('Communication versus inducing-point placement (M=%d, DAC rounds=%d)', ...
    m_value,dac_rounds));
exportgraphics(fig,fullfile(result_folder, ...
    sprintf('trigger_location_comparison_M%d_R%d.png',m_value,dac_rounds)), ...
    'Resolution',220);
savefig(fig,fullfile(result_folder, ...
    sprintf('trigger_location_comparison_M%d_R%d.fig',m_value,dac_rounds)));
end

function plot_one_trigger_set(rate_ax,activation_ax,data,label)
t_step = data.t_set(2)-data.t_set(1);
window_count = ceil(data.t_set(end));
rate = nan(1,window_count);
activation = nan(1,window_count);
trigger_data = double(data.dac_trigger_count_per_agent_point_set);
agent_active = squeeze(any(trigger_data > 0,2));
for window_i = 1:window_count
    mask = data.t_set(1:end-1) >= window_i-1 & ...
        data.t_set(1:end-1) < window_i;
    duration = nnz(mask)*t_step;
    rate(window_i) = sum(trigger_data(:,:,mask),'all') / ...
        (size(trigger_data,1)*size(trigger_data,2)*duration);
    activation(window_i) = mean(agent_active(:,mask),'all');
end
window_centers = (0:window_count-1)+0.5;
plot(rate_ax,window_centers,rate,'LineWidth',1.8,'DisplayName',label);
plot(activation_ax,window_centers,activation,'LineWidth',1.8,'DisplayName',label);
end

function plot_trajectory_dynamics(data,result_folder,m_value,dac_rounds)
valid = isfinite(data.prediction_error_norm_vector);
t = data.t_set(valid);
agent_colors = lines(size(data.f_true_all_set,2));
for output_i = 1:size(data.f_true_all_set,1)
    fig = figure('Color','w','Position',[120 80 1250 720]);
    layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
    methods = {data.f_true_all_set,'True dynamics';data.f_direct_all_set,'Direct PoE'; ...
        data.f_ip_cen_all_set,'IP-CEN';data.f_hat_all_set,'IP-DAC'};
    axes_set = gobjects(4,1);
    for method_i = 1:4
        axes_set(method_i) = nexttile(layout); hold(axes_set(method_i),'on');
        values = methods{method_i,1};
        for agent_i = 1:size(values,2)
            plot(axes_set(method_i),t,squeeze(values(output_i,agent_i,valid)), ...
                'Color',agent_colors(agent_i,:),'LineWidth',1.15, ...
                'DisplayName',sprintf('Agent %d',agent_i));
        end
        title(axes_set(method_i),methods{method_i,2});
        ylabel(axes_set(method_i),sprintf('f_%d',output_i));
        grid(axes_set(method_i),'on'); box(axes_set(method_i),'on');
    end
    linkaxes(axes_set,'xy'); xlabel(layout,'Time (s)');
    title(layout,sprintf('Trajectory-box IP: unknown dynamics output %d',output_i));
    legend(axes_set(2),'Location','best');
    exportgraphics(fig,fullfile(result_folder, ...
        sprintf('trajectory_dynamics_output%d_M%d_R%d.png',output_i,m_value,dac_rounds)), ...
        'Resolution',220);
    savefig(fig,fullfile(result_folder, ...
        sprintf('trajectory_dynamics_output%d_M%d_R%d.fig',output_i,m_value,dac_rounds)));
end
end

function restore_environment(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
