function output_file = plot_tp_ac_round_convergence( ...
    result_file,output_folder,target_time)
%PLOT_TP_AC_ROUND_CONVERGENCE Show what the ten TP-AC rounds accomplish.

if nargin < 3 || isempty(target_time), target_time = 10; end
if nargin < 2 || isempty(output_folder)
    output_folder = fileparts(result_file);
end
if ~isfolder(output_folder), mkdir(output_folder); end

d = load(result_file,'t_set','tp_ac_round_error_set', ...
    'TPACFixedRounds','TPQueryUpdateInterval','CurrentMode');
if ~isfield(d,'tp_ac_round_error_set')
    error('Result was generated before TP AC-round diagnostics were added.');
end

update_columns = find(isfinite(d.tp_ac_round_error_set(1,:)));
if isempty(update_columns)
    error('No TP-AC reference updates were recorded in %s.',result_file);
end
[~,nearest_idx] = min(abs(d.t_set(update_columns)-target_time));
time_idx = update_columns(nearest_idx);
window_start = d.t_set(time_idx);
window_end = window_start+d.TPQueryUpdateInterval;

round_error = d.tp_ac_round_error_set(:,time_idx);
valid = isfinite(round_error);
rounds = (0:d.TPACFixedRounds).';
rounds = rounds(valid);
round_error = round_error(valid);
reduction = 100*(1-round_error(end)/max(round_error(1),eps));

fig = figure('Visible','off','Color','w','Position',[100 100 760 520]);
ax = axes(fig);
plot(ax,rounds,round_error,'-o','Color',[0.12 0.12 0.12], ...
    'MarkerFaceColor','w','MarkerEdgeColor',[0.12 0.12 0.12], ...
    'LineWidth',1.6,'MarkerSize',5);
grid(ax,'on'); box(ax,'on');
xlim(ax,[0 d.TPACFixedRounds]);
xticks(ax,0:d.TPACFixedRounds);
ylim(ax,[0 1.08*max(round_error)]);
xlabel(ax,'AC round');
ylabel(ax,'RMS error to current average P');
title(ax,sprintf(['Selected window %.2f--%.2f s: ' ...
    'what ten AC rounds accomplish'],window_start,window_end), ...
    'FontWeight','normal');

annotation_text = sprintf('%.3g \\rightarrow %.3g\n(%.1f%% reduction)', ...
    round_error(1),round_error(end),reduction);
text(ax,d.TPACFixedRounds-0.25,0.94*max(round_error),annotation_text, ...
    'HorizontalAlignment','right','VerticalAlignment','top', ...
    'Color',[0.72 0.16 0.12],'FontSize',10);
set(ax,'FontName','Times New Roman','FontSize',11,'LineWidth',0.8);

[~,result_name] = fileparts(result_file);
output_file = fullfile(output_folder, ...
    sprintf('%s_ac_round_convergence_t%05.2f.png',result_name,window_start));
exportgraphics(fig,output_file,'Resolution',220);
savefig(fig,strrep(output_file,'.png','.fig'));
close(fig);
fprintf('Saved TP-AC round-convergence plot: %s\n',output_file);
end
