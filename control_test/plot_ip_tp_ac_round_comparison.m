function output_file = plot_ip_tp_ac_round_comparison( ...
    ip_result_file,tp_result_file,output_folder,target_time)
%PLOT_IP_TP_AC_ROUND_COMPARISON Compare fixed-R AC convergence fairly.
% Each curve is normalized by its round-0 RMS because IP and TP consensus
% states have different dimensions and raw numerical scales.

if nargin < 4 || isempty(target_time), target_time = 10; end
if nargin < 3 || isempty(output_folder)
    output_folder = fileparts(tp_result_file);
end
if ~isfolder(output_folder), mkdir(output_folder); end

ip = load(ip_result_file,'t_set','ac_tracking_error_history_set', ...
    'ACFixedIterations','ACForceFullBroadcast');
tp = load(tp_result_file,'t_set','tp_ac_round_error_set', ...
    'TPACFixedRounds');

ip_candidates = find(~cellfun(@isempty,ip.ac_tracking_error_history_set));
[~,j] = min(abs(ip.t_set(ip_candidates)-target_time));
ip_time_idx = ip_candidates(j);
ip_error = ip.ac_tracking_error_history_set{ip_time_idx}(:);

tp_candidates = find(isfinite(tp.tp_ac_round_error_set(1,:)));
[~,j] = min(abs(tp.t_set(tp_candidates)-target_time));
tp_time_idx = tp_candidates(j);
tp_error = tp.tp_ac_round_error_set(:,tp_time_idx);
tp_error = tp_error(isfinite(tp_error));

if numel(ip_error) ~= numel(tp_error)
    error('IP and TP AC histories must contain the same number of rounds.');
end
rounds = (0:numel(ip_error)-1).';
ip_relative = ip_error/max(ip_error(1),eps);
tp_relative = tp_error/max(tp_error(1),eps);
ip_reduction = 100*(1-ip_relative(end));
tp_reduction = 100*(1-tp_relative(end));
if isfield(ip,'ACForceFullBroadcast') && ip.ACForceFullBroadcast
    ip_label = 'IP-AC (full communication)';
else
    ip_label = 'IP-AC (PETC)';
end

fig = figure('Visible','off','Color','w','Position',[100 100 780 520]);
ax = axes(fig);
semilogy(ax,rounds,ip_relative,'-o','Color',[0.12 0.12 0.12], ...
    'MarkerFaceColor','w','LineWidth',1.6,'MarkerSize',5, ...
    'DisplayName',ip_label);
hold(ax,'on');
semilogy(ax,rounds,tp_relative,'--s','Color',[0.45 0.45 0.45], ...
    'MarkerFaceColor','w','LineWidth',1.6,'MarkerSize',5, ...
    'DisplayName','TP-AC (full communication)');
grid(ax,'on'); box(ax,'on');
xlim(ax,[0 max(rounds)]);
xticks(ax,rounds);
xlabel(ax,'AC round');
ylabel(ax,'Relative RMS consensus error');
title(ax,sprintf('IP-AC vs TP-AC at t = %.2f s (R = %d)', ...
    tp.t_set(tp_time_idx),numel(rounds)-1),'FontWeight','normal');
legend(ax,'Location','southwest','Box','off');
text(ax,0.97*rounds(end),0.30, ...
    sprintf('IP: %.1f%% reduction',ip_reduction), ...
    'HorizontalAlignment','right','VerticalAlignment','bottom', ...
    'Color',[0.12 0.12 0.12]);
text(ax,0.97*rounds(end),0.20, ...
    sprintf('TP: %.2f%% reduction',tp_reduction), ...
    'HorizontalAlignment','right','VerticalAlignment','bottom', ...
    'Color',[0.35 0.35 0.35]);
set(ax,'FontName','Times New Roman','FontSize',11,'LineWidth',0.8);

output_file = fullfile(output_folder,sprintf( ...
    'ip_ac_vs_tp_ac_round_convergence_t%05.2f.png', ...
    tp.t_set(tp_time_idx)));
exportgraphics(fig,output_file,'Resolution',220);
savefig(fig,strrep(output_file,'.png','.fig'));
close(fig);
fprintf('Saved IP/TP AC comparison: %s\n',output_file);
end
