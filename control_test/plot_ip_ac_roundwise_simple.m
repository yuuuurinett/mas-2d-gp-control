function output_file = plot_ip_ac_roundwise_simple(result_file,agent_i,window_start)
%PLOT_IP_AC_ROUNDWISE_SIMPLE Sequential 0.01-s detector trace.
%
% Top: all ten AC checks inside every 0.1-s physical window are unfolded
% in chronological order.  This gives one LHS curve and one RHS curve.
% Bottom: the ten checks in one selected physical window.

if nargin < 2 || isempty(agent_i), agent_i = 4; end
if nargin < 3 || isempty(window_start), window_start = 10; end

d = load(result_file,'t_set','ac_iteration_count_set', ...
    'ac_inner_trigger_instance_set','ac_trigger_measure_set', ...
    'ac_trigger_threshold_set','ACFixedIterations', ...
    'ACPeriodicSigma','NumInducingPoints');

rounds = d.ACFixedIterations;
update_idx = find(d.ac_iteration_count_set > 0);
physical_time = reshape(d.t_set(update_idx),1,[]);
lhs = double(squeeze(d.ac_trigger_measure_set( ...
    agent_i,1:rounds,update_idx)));
rhs = double(squeeze(d.ac_trigger_threshold_set( ...
    agent_i,1:rounds,update_idx)));
events = logical(squeeze(d.ac_inner_trigger_instance_set( ...
    agent_i,1:rounds,update_idx)));

if size(lhs,1) ~= rounds
    lhs = reshape(lhs,rounds,[]);
    rhs = reshape(rhs,rounds,[]);
    events = reshape(events,rounds,[]);
end

lhs_color = [0.24 0.18 0.48];
rhs_color = [0.82 0.43 0.06];
event_color = [0.72 0.08 0.12];

fig = figure('Visible','off','Color','w','Position',[60 60 1900 950]);
layout = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(layout); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
physical_window = median(diff(physical_time),'omitnan');
communication_step = physical_window/rounds;
check_time = repmat(physical_time,rounds,1) + ...
    repmat((0:rounds-1)'*communication_step,1,numel(physical_time));
check_time = check_time(:);
lhs_trace = lhs(:);
rhs_trace = rhs(:);
event_trace = events(:);
plot(ax1,check_time,lhs_trace,'-','Color',lhs_color, ...
    'LineWidth',0.85,'DisplayName','LHS');
plot(ax1,check_time,rhs_trace,'--','Color',rhs_color, ...
    'LineWidth',0.85,'DisplayName','RHS');
plot(ax1,check_time(event_trace),lhs_trace(event_trace),'*', ...
    'Color',event_color,'MarkerSize',3.2,'LineWidth',0.5, ...
    'DisplayName','Actual broadcast');
xlim(ax1,[check_time(1),check_time(end)]);
xlabel(ax1,sprintf('Communication time (s), check interval = %.4f s', ...
    communication_step));
ylabel(ax1,'Detector value');
title(ax1,sprintf(['Sequential detector trace; %d checks and %d ' ...
    'broadcasts'],numel(events),nnz(events)),'FontWeight','normal');
legend(ax1,'Location','best','Box','off');

[~,window_i] = min(abs(physical_time-window_start));
ax2 = nexttile(layout); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
round_axis = 1:rounds;
plot(ax2,round_axis,lhs(:,window_i),'-o','Color',lhs_color, ...
    'LineWidth',1.4,'MarkerSize',4,'DisplayName','LHS');
plot(ax2,round_axis,rhs(:,window_i),'--o','Color',rhs_color, ...
    'LineWidth',1.4,'MarkerSize',4,'DisplayName','RHS');
trigger_window = events(:,window_i);
plot(ax2,round_axis(trigger_window),lhs(trigger_window,window_i),'*', ...
    'Color',event_color,'MarkerSize',10,'LineWidth',1.2, ...
    'DisplayName','Broadcast: LHS >= RHS');
xlim(ax2,[1 rounds]);
xticks(ax2,round_axis);
xlabel(ax2,'AC round, ell');
ylabel(ax2,'Detector value');
title(ax2,sprintf('Physical window starting at t = %.2f s', ...
    physical_time(window_i)),'FontWeight','normal');
legend(ax2,'Location','best','Box','off');

title(layout,sprintf(['Agent %d periodic ET-AC detector ' ...
    '(M=%d, R=%d, sigma=%.1f)'],agent_i,d.NumInducingPoints, ...
    rounds,d.ACPeriodicSigma),'FontWeight','bold');

[result_folder,result_stem] = fileparts(result_file);
output_file = fullfile(result_folder,sprintf( ...
    '%s_sequential_detector_agent%d_window_%0.1fs.png', ...
    result_stem,agent_i,physical_time(window_i)));
exportgraphics(fig,output_file,'Resolution',240);
exportgraphics(fig,strrep(output_file,'.png','.pdf'),'ContentType','vector');
savefig(fig,strrep(output_file,'.png','.fig'));
close(fig);

fprintf('Saved round-wise AC detector figure: %s\n',output_file);
end
