function summary = plot_prediction_error_bound(result_file,output_file)
%PLOT_PREDICTION_ERROR_BOUND Compare per-agent true error and GP bound.
%   SUMMARY = PLOT_PREDICTION_ERROR_BOUND(RESULT_FILE) loads one saved
%   inducing-point online simulation and plots
%       e_i(t) = ||f(x_i(t))-mu_com,i(x_i(t))||_2
%   against the three-layer bound assembled with beta=3 and gamma=0.01.

arguments
    result_file (1,1) string
    output_file (1,1) string = ""
end

required = {'t_set','prediction_error_per_agent_set', ...
    'prediction_error_bound_set','PredictionBoundBeta', ...
    'PredictionBoundGamma'};
S = load(result_file,required{:});
for idx = 1:numel(required)
    if ~isfield(S,required{idx})
        error('Result file is missing required variable "%s".',required{idx});
    end
end

t = S.t_set(:).';
error_set = S.prediction_error_per_agent_set;
bound_set = S.prediction_error_bound_set;
if ~isequal(size(error_set),size(bound_set)) || ...
        size(error_set,2) ~= numel(t)
    error('Saved error/bound arrays are inconsistent with t_set.');
end

valid = isfinite(error_set) & isfinite(bound_set);
covered = valid & error_set <= bound_set;
agent_count = size(error_set,1);
coverage_per_agent = nan(agent_count,1);
max_error = nan(agent_count,1);
max_bound = nan(agent_count,1);
min_margin = nan(agent_count,1);
for agent_i = 1:agent_count
    agent_valid = valid(agent_i,:);
    if any(agent_valid)
        coverage_per_agent(agent_i) = mean(covered(agent_i,agent_valid));
        max_error(agent_i) = max(error_set(agent_i,agent_valid));
        max_bound(agent_i) = max(bound_set(agent_i,agent_valid));
        min_margin(agent_i) = min( ...
            bound_set(agent_i,agent_valid)-error_set(agent_i,agent_valid));
    end
end
summary = table((1:agent_count).',coverage_per_agent,max_error,max_bound,min_margin, ...
    'VariableNames',{'Agent','Coverage','MaxError','MaxBound','MinMargin'});
disp(summary);
fprintf('Overall empirical coverage: %.2f%%\n',100*mean(covered(valid)));

fig = figure('Color','w','Position',[100 100 1150 720]);
layout = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
title(layout,sprintf('Prediction error bound (\\beta=%.1f, \\gamma=%.2f)', ...
    S.PredictionBoundBeta,S.PredictionBoundGamma),'Interpreter','tex');

colors = lines(agent_count);
ax1 = nexttile(layout,1);
hold(ax1,'on');
legend_handles = gobjects(2*agent_count,1);
legend_text = strings(2*agent_count,1);
for agent_i = 1:agent_count
    legend_handles(2*agent_i-1) = plot(ax1,t,error_set(agent_i,:),'-', ...
        'Color',colors(agent_i,:),'LineWidth',1.2);
    legend_handles(2*agent_i) = plot(ax1,t,bound_set(agent_i,:),'--', ...
        'Color',colors(agent_i,:),'LineWidth',1.2);
    legend_text(2*agent_i-1) = sprintf('Agent %d error',agent_i);
    legend_text(2*agent_i) = sprintf('Agent %d bound',agent_i);
end
grid(ax1,'on');
xlabel(ax1,'Time (s)');
ylabel(ax1,'Prediction error / bound');
legend(ax1,legend_handles,legend_text,'Location','eastoutside', ...
    'NumColumns',2);

ax2 = nexttile(layout,2);
hold(ax2,'on');
for agent_i = 1:agent_count
    plot(ax2,t,bound_set(agent_i,:)-error_set(agent_i,:), ...
        'Color',colors(agent_i,:),'LineWidth',1.1);
end
yline(ax2,0,'k--','LineWidth',1.2);
grid(ax2,'on');
xlabel(ax2,'Time (s)');
ylabel(ax2,'Bound - true error');
title(ax2,'Nonnegative margin means the bound is satisfied');
linkaxes([ax1 ax2],'x');

if strlength(output_file) == 0
    [folder,name] = fileparts(result_file);
    output_file = fullfile(folder,name + "_prediction_error_bound.png");
end
exportgraphics(fig,output_file,'Resolution',300);
fprintf('Saved figure: %s\n',output_file);
end
