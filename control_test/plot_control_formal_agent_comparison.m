function plot_control_formal_agent_comparison(ResultFolder)
% Formal 0.1 s online-update comparison. Each agent has one color.

if nargin < 1 || isempty(ResultFolder)
    ProjectRoot = fileparts(fileparts(mfilename('fullpath')));
    ResultFolder = fullfile(ProjectRoot,'result','control_formal_0p1');
end

C = load(fullfile(ResultFolder,'centralized_poe_online_0p1.mat'));
N = load(fullfile(ResultFolder,'neighbor_poe_online_0p1.mat'));
L = load(fullfile(ResultFolder,'local_online_0p1.mat'));
E = load(fullfile(ResultFolder,'exact_reference.mat'));
sets = {C,N};
names = {'Centralized PoE','Neighbor PoE'};
colors = lines(size(C.f_true_all_set,2));

%% 1. Per-agent tracking error
tracking_sets = {C,N,L,E};
tracking_names = {'Centralized PoE','Neighbor PoE','Local GP','Exact dynamics'};
fig = figure('Color','w','Position',[40,35,1400,900]);
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for method_i = 1:numel(tracking_sets)
    S = tracking_sets{method_i};
    ax = nexttile(layout); hold(ax,'on');
    tracking = per_agent_tracking(S.vartheta_all_set);
    for agent_i = 1:size(tracking,1)
        semilogy(ax,S.t_set,max(tracking(agent_i,:),eps),'-', ...
            'Color',colors(agent_i,:),'LineWidth',1.4, ...
            'DisplayName',sprintf('Agent %d',agent_i));
    end
    format_axis(ax,'Time (s)','Tracking error norm',tracking_names{method_i});
    legend(ax,'Location','best');
end
title(layout,'Per-agent tracking error (online update every 0.1 s)');
export_clean(fig);
save_both(fig,ResultFolder,'formal_tracking_error_per_agent');

%% 2. Predicted versus true dynamics, vector norm over both outputs
fig = figure('Color','w','Position',[40,35,1400,900]);
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for method_i = 1:2
    S = sets{method_i};
    prediction_ax = nexttile(layout,2*method_i-1); hold(prediction_ax,'on');
    true_ax = nexttile(layout,2*method_i); hold(true_ax,'on');
    all_norms = [];
    for agent_i = 1:size(S.f_true_all_set,2)
        true_norm = squeeze(vecnorm(S.f_true_all_set(:,agent_i,:),2,1));
        pred_norm = squeeze(vecnorm(S.f_hat_all_set(:,agent_i,:),2,1));
        all_norms = [all_norms; true_norm(:); pred_norm(:)]; %#ok<AGROW>
        plot(prediction_ax,S.t_set,pred_norm,'-', ...
            'Color',colors(agent_i,:),'LineWidth',1.2, ...
            'DisplayName',sprintf('Agent %d',agent_i));
        plot(true_ax,S.t_set,true_norm,'-', ...
            'Color',colors(agent_i,:),'LineWidth',1.2, ...
            'DisplayName',sprintf('Agent %d',agent_i));
    end
    common_ylim = [min(all_norms),max(all_norms)];
    common_ylim = common_ylim + [-1,1]*0.03*max(diff(common_ylim),eps);
    ylim(prediction_ax,common_ylim); ylim(true_ax,common_ylim);
    format_axis(prediction_ax,'Time (s)','Dynamics vector norm', ...
        [names{method_i},': prediction']);
    format_axis(true_ax,'Time (s)','Dynamics vector norm', ...
        [names{method_i},': true dynamics']);
    legend(prediction_ax,'Location','best');
end
title(layout,'Prediction and true unknown dynamics (two outputs combined by vector norm)');
export_clean(fig);
save_both(fig,ResultFolder,'formal_prediction_vs_true_per_agent');

%% 3. Prediction error versus tracking error on a shared time axis
fig = figure('Color','w','Position',[50,35,1400,850]);
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for method_i = 1:2
    S = sets{method_i};
    tracking = per_agent_tracking(S.vartheta_all_set);
    ax_prediction = nexttile(layout,method_i); hold(ax_prediction,'on');
    ax_tracking = nexttile(layout,method_i+2); hold(ax_tracking,'on');
    for agent_i = 1:size(tracking,1)
        prediction_error = squeeze(vecnorm( ...
            S.f_hat_all_set(:,agent_i,:)-S.f_true_all_set(:,agent_i,:),2,1));
        semilogy(ax_prediction,S.t_set,max(prediction_error,eps),'-', ...
            'Color',colors(agent_i,:),'LineWidth',1.2, ...
            'DisplayName',sprintf('Agent %d',agent_i));
        semilogy(ax_tracking,S.t_set,max(tracking(agent_i,:),eps),'-', ...
            'Color',colors(agent_i,:),'LineWidth',1.2, ...
            'DisplayName',sprintf('Agent %d',agent_i));
    end
    format_axis(ax_prediction,'Time (s)','Prediction-error norm', ...
        [names{method_i},': prediction error']);
    format_axis(ax_tracking,'Time (s)','Tracking-error norm', ...
        [names{method_i},': tracking error']);
    legend(ax_prediction,'Location','best');
end
title(layout,'Prediction error versus tracking error on the same time axis');
export_clean(fig);
save_both(fig,ResultFolder,'formal_prediction_error_vs_tracking_error_per_agent');
end

function tracking = per_agent_tracking(vartheta)
AgentQuantity = 6;
x_dim = size(vartheta,1)/AgentQuantity;
tracking = zeros(AgentQuantity,size(vartheta,2));
for agent_i = 1:AgentQuantity
    idx = (agent_i-1)*x_dim+(1:x_dim);
    tracking(agent_i,:) = vecnorm(vartheta(idx,:),2,1);
end
end

function format_axis(ax,x_text,y_text,title_text)
grid(ax,'on'); box(ax,'on');
xlabel(ax,x_text); ylabel(ax,y_text); title(ax,title_text);
end

function export_clean(fig)
axes_set = findall(fig,'Type','axes');
for ax_i = 1:numel(axes_set)
    disableDefaultInteractivity(axes_set(ax_i));
    set(axes_set(ax_i),'Toolbar',[]);
end
end

function save_both(fig,folder,stem)
exportgraphics(fig,fullfile(folder,[stem,'.png']),'Resolution',220);
savefig(fig,fullfile(folder,[stem,'.fig']));
end
