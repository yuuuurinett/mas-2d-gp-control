function output_file = plot_ip_ac_agent_level_two_lines( ...
    result_file,agent_i,end_time)
%PLOT_IP_AC_AGENT_LEVEL_TWO_LINES Wide, non-averaged PETC detector plot.
% Every stored communication check is plotted.  The two displayed signals
% are exactly the squared agent-level detector sides
%   LHS_i = ||Xi_i-XiHat_i||_F^2/(pM),
%   RHS_i = sigma_i^2 sum_j a_ij||XiHat_i-XiHat_j||_F^2/(pM).

if nargin < 2 || isempty(agent_i), agent_i = 4; end
if nargin < 3 || isempty(end_time), end_time = 10; end

d = load(result_file,'t_set','ac_iteration_count_set', ...
    'ac_trigger_measure_set','ac_trigger_threshold_set', ...
    'ac_inner_trigger_instance_set','ACFixedIterations', ...
    'projection_update_set');

available_variables = who('-file',result_file);
detector_is_squared = false;
if any(strcmp(available_variables,'ACDetectorStoresSquared'))
    metadata = load(result_file,'ACDetectorStoresSquared');
    detector_is_squared = logical(metadata.ACDetectorStoresSquared);
end

update_idx = find(d.ac_iteration_count_set > 0);
window_time = d.t_set(update_idx);
is_physical_time_ac = any(strcmp(available_variables, ...
    'ACCommunicationInterval'));
if is_physical_time_ac
    timing = load(result_file,'ACCommunicationInterval', ...
        'ACChecksPerPhysicalStep');
    delta_c = timing.ACCommunicationInterval;
    checks_per_step = timing.ACChecksPerPhysicalStep;
else
    checks_per_step = d.ACFixedIterations;
    delta_c = 0.1/checks_per_step;
end

lhs = double(squeeze(d.ac_trigger_measure_set( ...
    agent_i,1:checks_per_step,update_idx)));
rhs = double(squeeze(d.ac_trigger_threshold_set( ...
    agent_i,1:checks_per_step,update_idx)));
events = logical(squeeze(d.ac_inner_trigger_instance_set( ...
    agent_i,1:checks_per_step,update_idx)));
if checks_per_step == 1
    lhs = reshape(lhs,1,[]);
    rhs = reshape(rhs,1,[]);
    events = reshape(events,1,[]);
end
if ~detector_is_squared
    % Legacy cache: sqrt(LHS) and sqrt(RHS) were stored. Squaring restores
    % the paper's literal detector without changing any trigger decision.
    lhs = lhs.^2;
    rhs = rhs.^2;
end

time = repmat(window_time,checks_per_step,1) + ...
    repmat((0:checks_per_step-1)'*delta_c,1,numel(window_time));

% Interpolate only genuine sign-bracketed crossings.  No detector samples
% or broadcasts are removed from the two main curves.
cross_time = [];
cross_value = [];
time_sequence = time(:);
lhs_sequence = lhs(:);
rhs_sequence = rhs(:);
event_sequence = events(:);
reset_sequence = repelem(logical(d.projection_update_set(update_idx)), ...
    checks_per_step).';
g_sequence = lhs_sequence-rhs_sequence;
for sample_i = 2:numel(time_sequence)
    if ~reset_sequence(sample_i) && event_sequence(sample_i) && ...
            isfinite(g_sequence(sample_i-1)) && ...
            isfinite(g_sequence(sample_i)) && ...
            g_sequence(sample_i-1) < 0 && g_sequence(sample_i) >= 0
        fraction = -g_sequence(sample_i-1) / ...
            (g_sequence(sample_i)-g_sequence(sample_i-1));
        cross_time(end+1,1) = time_sequence(sample_i-1) + ... %#ok<AGROW>
            fraction*(time_sequence(sample_i)-time_sequence(sample_i-1));
        lhs_cross = lhs_sequence(sample_i-1)+fraction* ...
            (lhs_sequence(sample_i)-lhs_sequence(sample_i-1));
        rhs_cross = rhs_sequence(sample_i-1)+fraction* ...
            (rhs_sequence(sample_i)-rhs_sequence(sample_i-1));
        cross_value(end+1,1) = 0.5*(lhs_cross+rhs_cross); %#ok<AGROW>
    end
end

time_vector = time(:);
lhs_vector = lhs(:);
rhs_vector = rhs(:);
display_mask = time_vector <= end_time;
cross_mask = cross_time <= end_time;

fig = figure('Visible','off','Color','w', ...
    'Position',[40 80 6000 900]);
ax = axes(fig,'Position',[0.055 0.14 0.925 0.78]);
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax,time_vector(display_mask),lhs_vector(display_mask),'-', ...
    'Color',[0.16 0.34 0.66],'LineWidth',0.75, ...
    'DisplayName','LHS_i');
plot(ax,time_vector(display_mask),rhs_vector(display_mask),'-', ...
    'Color',[0.85 0.33 0.10],'LineWidth',0.75, ...
    'DisplayName','RHS_i');
plot(ax,cross_time(cross_mask),cross_value(cross_mask),'kx', ...
    'MarkerSize',4.5,'LineWidth',1.0, ...
    'DisplayName','LHS_i = RHS_i crossing');
xlim(ax,[0 end_time]);
xlabel(ax,'Communication time (s)');
ylabel(ax,'Squared detector value (normalized by pM)');
title(ax,sprintf(['Agent %d periodic ET-AC, \\Delta_c=%.4f s: ' ...
    'every check shown, no averaging'],agent_i,delta_c), ...
    'FontWeight','normal');
legend(ax,'Location','northoutside','Orientation','horizontal', ...
    'Box','off');
set(ax,'FontSize',12,'Layer','top');

output_folder = fileparts(result_file);
output_file = fullfile(output_folder,sprintf( ...
    'agent%d_two_detector_lines_T%d.png',agent_i,round(end_time)));
saveas(fig,output_file);
savefig(fig,strrep(output_file,'.png','.fig'));
close(fig);
end
