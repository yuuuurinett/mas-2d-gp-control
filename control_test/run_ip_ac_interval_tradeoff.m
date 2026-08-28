function Summary = run_ip_ac_interval_tradeoff(do_run)
%RUN_IP_AC_INTERVAL_TRADEOFF Periodic ET-AC time-resolution comparison.
% P is refreshed every 0.1 s. Within each frozen-P window, AC checks the
% detector at DeltaC and uses h=DeltaC*KappaP, so every case approximates
% the same continuous consensus gain KappaP.

if nargin < 1 || isempty(do_run), do_run = true; end
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result', ...
    'ip_ac_interval_tradeoff_u0p6_T20');
if ~isfolder(result_folder), mkdir(result_folder); end

delta_c = [0.01;0.005;0.0025];
kappa_p = 20;
rounds = round(0.1./delta_c);
consensus_step = delta_c*kappa_p;
case_count = numel(delta_c);
result_files = strings(case_count,1);

point_file = fullfile(repo_root,'result', ...
    'inducing_points_position_grid100_velocity6_full_domain.mat');
env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_IP_INDUCING_POINT_FILE','CONTROL_IP_POSITION_LENGTH_SCALE', ...
    'CONTROL_IP_VELOCITY_LENGTH_SCALE','CONTROL_AC_FIXED_ITERATIONS', ...
    'CONTROL_AC_PERIODIC_SIGMA','CONTROL_AC_CONSENSUS_STEP', ...
    'CONTROL_AC_TRIGGER_DIAGNOSTICS'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment_local( ...
    env_names,old_values)); %#ok<NASGU>
setenv('CONTROL_SIM_END_TIME','20');
setenv('CONTROL_IP_NUM_INDUCING_POINTS','600');
setenv('CONTROL_IP_INDUCING_POINT_FILE',point_file);
setenv('CONTROL_IP_POSITION_LENGTH_SCALE','0.5');
setenv('CONTROL_IP_VELOCITY_LENGTH_SCALE','1.5');
setenv('CONTROL_AC_PERIODIC_SIGMA','0.3');
setenv('CONTROL_AC_TRIGGER_DIAGNOSTICS','1');

for case_i = 1:case_count
    setenv('CONTROL_AC_FIXED_ITERATIONS',num2str(rounds(case_i)));
    setenv('CONTROL_AC_CONSENSUS_STEP', ...
        num2str(consensus_step(case_i),'%.12g'));
    delta_tag = strrep(num2str(delta_c(case_i),'%.4f'),'.','p');
    name = sprintf('ip_ac_dc_%s_R%d_h_%0.3f_seed42', ...
        delta_tag,rounds(case_i),consensus_step(case_i));
    name = strrep(name,'.','p');
    result_files(case_i) = fullfile(result_folder,[name '.mat']);
    if do_run && ~isfile(result_files(case_i))
        run_simulation_inducing_point('poe_ac',result_folder,name, ...
            false,0.6,0.1,[],42);
    end
end

MaxTrackingAfter10 = nan(case_count,1);
MeanTrackingAfter10 = nan(case_count,1);
MeanPredictionAfter10 = nan(case_count,1);
BroadcastsPerAgent = nan(case_count,1);
MedianTriggerOvershoot = nan(case_count,1);
MedianRelativeOvershoot = nan(case_count,1);
MeanTerminalACError = nan(case_count,1);

for case_i = 1:case_count
    d = load(result_files(case_i));
    tracking = sqrt(sum(d.vartheta_all_set.^2,1));
    stable = d.t_set >= 10;
    MaxTrackingAfter10(case_i) = max(tracking(stable),[],'omitnan');
    MeanTrackingAfter10(case_i) = mean(tracking(stable),'omitnan');
    MeanPredictionAfter10(case_i) = mean( ...
        d.prediction_error_norm_vector(stable),'omitnan');
    BroadcastsPerAgent(case_i) = d.ac_broadcasts_per_agent;

    update_idx = find(d.ac_iteration_count_set > 0);
    lhs = double(d.ac_trigger_measure_set(:,1:rounds(case_i),update_idx));
    rhs = double(d.ac_trigger_threshold_set(:,1:rounds(case_i),update_idx));
    % Older cached runs stored the square-root-equivalent detector. Convert
    % them to the paper's squared LHS/RHS before reporting or plotting.
    if ~isfield(d,'ACDetectorStoresSquared') || ...
            ~logical(d.ACDetectorStoresSquared)
        lhs = lhs.^2;
        rhs = rhs.^2;
    end
    events = logical(d.ac_inner_trigger_instance_set( ...
        :,1:rounds(case_i),update_idx));
    gap = lhs(events)-rhs(events);
    event_rhs = rhs(events);
    valid = isfinite(gap) & isfinite(event_rhs) & event_rhs > 1e-10;
    MedianTriggerOvershoot(case_i) = median(gap(valid),'omitnan');
    MedianRelativeOvershoot(case_i) = median( ...
        gap(valid)./event_rhs(valid),'omitnan');

    histories = d.ac_tracking_error_history_set(update_idx);
    terminal_error = cellfun(@(x)x(end),histories);
    MeanTerminalACError(case_i) = mean(terminal_error,'omitnan');

    plot_interval_detector_local(result_files(case_i),delta_c(case_i), ...
        4,10,0.5,result_folder);
end

DeltaC = delta_c;
KappaP = kappa_p*ones(case_count,1);
RoundsPerPointUpdate = rounds;
ConsensusStep = consensus_step;
Summary = table(DeltaC,KappaP,RoundsPerPointUpdate,ConsensusStep, ...
    MaxTrackingAfter10,MeanTrackingAfter10,MeanPredictionAfter10, ...
    BroadcastsPerAgent,MedianTriggerOvershoot,MedianRelativeOvershoot, ...
    MeanTerminalACError,result_files);
writetable(Summary,fullfile(result_folder,'interval_tradeoff_summary.csv'));
save(fullfile(result_folder,'interval_tradeoff_summary.mat'),'Summary');
disp(Summary(:,1:11));
end

function plot_interval_detector_local(result_file,delta_c,agent_i, ...
    zoom_start,zoom_span,output_folder)
d = load(result_file,'t_set','ac_iteration_count_set', ...
    'ac_trigger_measure_set','ac_trigger_threshold_set', ...
    'ac_inner_trigger_instance_set','ACFixedIterations');
available_variables = who('-file',result_file);
detector_is_squared = false;
if any(strcmp(available_variables,'ACDetectorStoresSquared'))
    detector_metadata = load(result_file,'ACDetectorStoresSquared');
    detector_is_squared = logical( ...
        detector_metadata.ACDetectorStoresSquared);
end
rounds = d.ACFixedIterations;
update_idx = find(d.ac_iteration_count_set > 0);
window_time = d.t_set(update_idx);
lhs_matrix = double(squeeze(d.ac_trigger_measure_set( ...
    agent_i,1:rounds,update_idx)));
rhs_matrix = double(squeeze(d.ac_trigger_threshold_set( ...
    agent_i,1:rounds,update_idx)));
if ~detector_is_squared
    lhs_matrix = lhs_matrix.^2;
    rhs_matrix = rhs_matrix.^2;
end
event_matrix = logical(squeeze(d.ac_inner_trigger_instance_set( ...
    agent_i,1:rounds,update_idx)));
if rounds == 1
    lhs_matrix = reshape(lhs_matrix,1,[]);
    rhs_matrix = reshape(rhs_matrix,1,[]);
    event_matrix = reshape(event_matrix,1,[]);
end
time_matrix = repmat(window_time,rounds,1) + ...
    repmat((0:rounds-1)'*delta_c,1,numel(window_time));

time_all = time_matrix(:);
lhs_all = lhs_matrix(:);
rhs_all = rhs_matrix(:);
event_all = event_matrix(:);

% Reconstruct continuous crossing locations only when the sampled detector
% brackets g=LHS-RHS=0 inside the same frozen-P window. Actual statistics
% continue to use every sampled broadcast event.
cross_time = [];
cross_value = [];
for window_i = 1:numel(window_time)
    g = lhs_matrix(:,window_i)-rhs_matrix(:,window_i);
    for round_i = 2:rounds
        if event_matrix(round_i,window_i) && ...
                isfinite(g(round_i-1)) && isfinite(g(round_i)) && ...
                g(round_i-1) < 0 && g(round_i) >= 0
            fraction = -g(round_i-1)/(g(round_i)-g(round_i-1));
            t_cross = time_matrix(round_i-1,window_i) + ...
                fraction*delta_c;
            lhs_cross = lhs_matrix(round_i-1,window_i) + ...
                fraction*(lhs_matrix(round_i,window_i) ...
                -lhs_matrix(round_i-1,window_i));
            rhs_cross = rhs_matrix(round_i-1,window_i) + ...
                fraction*(rhs_matrix(round_i,window_i) ...
                -rhs_matrix(round_i-1,window_i));
            cross_time(end+1,1) = t_cross; %#ok<AGROW>
            cross_value(end+1,1) = 0.5*(lhs_cross+rhs_cross); %#ok<AGROW>
        end
    end
end

fig = figure('Visible','off','Color','w','Position',[80 80 1450 760]);
layout = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(layout); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
% Insert NaNs between GP windows.  This keeps every actual detector sample
% without drawing an artificial line from the final AC round of one window
% to the reset at the start of the next window.
time_segmented = [time_matrix;nan(1,numel(window_time))];
lhs_segmented = [lhs_matrix;nan(1,numel(window_time))];
rhs_segmented = [rhs_matrix;nan(1,numel(window_time))];
plot(ax1,time_segmented(:),lhs_segmented(:),'-', ...
    'Color',[0.18 0.36 0.68],'LineWidth',0.65, ...
    'DisplayName','LHS: local broadcast error');
plot(ax1,time_segmented(:),rhs_segmented(:),'-', ...
    'Color',[0.82 0.36 0.12],'LineWidth',0.65, ...
    'DisplayName','RHS: neighbor threshold');
max_markers = 180;
marker_step = max(1,ceil(numel(cross_time)/max_markers));
display_idx = 1:marker_step:numel(cross_time);
plot(ax1,cross_time(display_idx),cross_value(display_idx),'x', ...
    'Color',[0.10 0.10 0.10],'MarkerSize',5,'LineWidth',1.0, ...
    'DisplayName','Interpolated LHS=RHS crossing');
xlim(ax1,[0,max(time_all)]);
xlabel(ax1,'Communication time (s)');
ylabel(ax1,'Squared detector value (normalized by pM)');
title(ax1,sprintf(['Agent %d, \\Delta_c=%.4f s: all actual checks, ' ...
    'no averaging (%d broadcasts; crossing markers shown 1/%d)'], ...
    agent_i,delta_c,nnz(event_all),marker_step),'FontWeight','normal');
legend(ax1,'Location','best','Box','off');

ax2 = nexttile(layout); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
zoom_end = zoom_start+zoom_span;
zoom_mask = time_all >= zoom_start & time_all <= zoom_end;
plot(ax2,time_all(zoom_mask),lhs_all(zoom_mask),'-o', ...
    'Color',[0.18 0.36 0.68],'MarkerSize',2.5,'LineWidth',1.1, ...
    'DisplayName','LHS');
plot(ax2,time_all(zoom_mask),rhs_all(zoom_mask),'-', ...
    'Color',[0.82 0.36 0.12],'LineWidth',1.1,'DisplayName','RHS');
cross_zoom = cross_time >= zoom_start & cross_time <= zoom_end;
plot(ax2,cross_time(cross_zoom),cross_value(cross_zoom),'kx', ...
    'MarkerSize',7,'LineWidth',1.3,'DisplayName','LHS=RHS crossing');
xlim(ax2,[zoom_start zoom_end]);
xlabel(ax2,'Communication time (s)');
ylabel(ax2,'Squared detector value (normalized by pM)');
title(ax2,sprintf('Local enlargement %.2f--%.2f s: all crossings shown', ...
    zoom_start,zoom_end),'FontWeight','normal');
legend(ax2,'Location','best','Box','off');

[~,base_name] = fileparts(result_file); %#ok<ASGLU>
delta_tag = strrep(sprintf('%.4f',delta_c),'.','p');
output_file = fullfile(char(output_folder), ...
    ['agent4_detector_agent_level_squared_dc_' delta_tag '.png']);
saveas(fig,output_file);
savefig(fig,strrep(output_file,'.png','.fig'));
close(fig);
end

function restore_environment_local(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
