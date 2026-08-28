function Summary = plot_ip_agent_m400_grid()
%PLOT_IP_AGENT_M400_GRID Summarize M=400 DAC rounds/epsilon grid.

repo_root = fileparts(fileparts(mfilename('fullpath')));
result_folder = fullfile(repo_root,'result','ip_agent_m400_grid');
round_values = [1 2 5];
epsilon_values = [0.04 0.08 0.1];
case_count = numel(round_values)*numel(epsilon_values);
data_set = cell(numel(round_values),numel(epsilon_values));

R = nan(case_count,1); Epsilon = nan(case_count,1);
MeanBroadcastRate = nan(case_count,1);
First5sBroadcastRate = nan(case_count,1);
Last5sBroadcastRate = nan(case_count,1);
LateRateSlope = nan(case_count,1);
MeanTrackingError = nan(case_count,1);
FinalTrackingError = nan(case_count,1);
MeanPredictionError = nan(case_count,1);
MeanDACTrackingError = nan(case_count,1);
LateDACTrackingError = nan(case_count,1);

rate_fig = figure('Color','w','Position',[80 40 1320 900]);
rate_layout = tiledlayout(rate_fig,3,3, ...
    'TileSpacing','compact','Padding','compact');
case_i = 0;
for row_i = 1:numel(round_values)
    for col_i = 1:numel(epsilon_values)
        case_i = case_i+1;
        round_count = round_values(row_i);
        epsilon_value = epsilon_values(col_i);
        epsilon_tag = strrep(sprintf('%.3g',epsilon_value),'.','p');
        result_file = fullfile(result_folder,sprintf( ...
            'poe_M400_R%d_eps_%s_T30.mat',round_count,epsilon_tag));
        data = load(result_file);
        data_set{row_i,col_i} = data;

        events = double(data.dac_inner_trigger_instance_set);
        event_t = data.t_set(1:size(events,3));
        duration = data.t_set(end)-data.t_set(1);
        agent_count = size(events,1);
        window_count = ceil(duration);
        rate = nan(1,window_count);
        for window_i = 1:window_count
            mask = event_t >= window_i-1 & event_t < window_i;
            rate(window_i) = sum(events(:,:,mask),'all')/agent_count;
        end
        centers = (0:window_count-1)+0.5;

        R(case_i) = round_count;
        Epsilon(case_i) = epsilon_value;
        MeanBroadcastRate(case_i) = mean(rate);
        First5sBroadcastRate(case_i) = mean(rate(1:5));
        Last5sBroadcastRate(case_i) = mean(rate(end-4:end));
        late_mask = centers >= 10;
        late_fit = polyfit(centers(late_mask),rate(late_mask),1);
        LateRateSlope(case_i) = late_fit(1);
        MeanTrackingError(case_i) = trapz( ...
            data.t_set,data.TrackingError_vector)/duration;
        FinalTrackingError(case_i) = data.TrackingError_vector(end);
        MeanPredictionError(case_i) = mean( ...
            data.prediction_error_norm_vector,'omitnan');
        MeanDACTrackingError(case_i) = mean( ...
            data.dac_tracking_error_set,'omitnan');
        late_dac_mask = event_t >= duration-5;
        LateDACTrackingError(case_i) = mean( ...
            data.dac_tracking_error_set(late_dac_mask),'omitnan');

        ax = nexttile(rate_layout); hold(ax,'on');
        plot(ax,centers,rate,'k-','LineWidth',1.25);
        yline(ax,Last5sBroadcastRate(case_i),'--', ...
            'Color',[0.85 0.33 0.1],'LineWidth',1.0);
        xlim(ax,[0 duration]); ylim(ax,[0 inf]);
        grid(ax,'on'); box(ax,'on');
        title(ax,sprintf('R=%d, \\epsilon=%.2f',round_count,epsilon_value));
        if col_i == 1, ylabel(ax,'Broadcasts/agent/s'); end
        if row_i == numel(round_values), xlabel(ax,'Time (s)'); end
    end
end
title(rate_layout,'M=400 IP-DAC packaged broadcasts (1-s windows)');
exportgraphics(rate_fig,fullfile(result_folder, ...
    'm400_broadcast_rate_grid.png'),'Resolution',230);
savefig(rate_fig,fullfile(result_folder,'m400_broadcast_rate_grid.fig'));

Summary = table(R,Epsilon,MeanBroadcastRate,First5sBroadcastRate, ...
    Last5sBroadcastRate,LateRateSlope,MeanTrackingError, ...
    FinalTrackingError,MeanPredictionError,MeanDACTrackingError, ...
    LateDACTrackingError);
writetable(Summary,fullfile(result_folder,'summary.csv'));
save(fullfile(result_folder,'summary.mat'),'Summary');

tradeoff_fig = figure('Color','w','Position',[100 100 1120 470]);
tradeoff_layout = tiledlayout(tradeoff_fig,1,2, ...
    'TileSpacing','compact','Padding','compact');
metrics = {MeanTrackingError,MeanPredictionError};
ylabels = {'Mean tracking error','Mean prediction error'};
for panel_i = 1:2
    ax = nexttile(tradeoff_layout); hold(ax,'on');
    for case_i = 1:height(Summary)
        scatter(ax,MeanBroadcastRate(case_i),metrics{panel_i}(case_i), ...
            55,R(case_i),'filled');
        text(ax,MeanBroadcastRate(case_i),metrics{panel_i}(case_i), ...
            sprintf('  R%d,%.2f',R(case_i),Epsilon(case_i)), ...
            'FontSize',8,'VerticalAlignment','middle');
    end
    xlabel(ax,'Mean broadcasts/agent/s'); ylabel(ax,ylabels{panel_i});
    grid(ax,'on'); box(ax,'on');
end
title(tradeoff_layout,'M=400 communication-accuracy tradeoff');
exportgraphics(tradeoff_fig,fullfile(result_folder, ...
    'm400_accuracy_tradeoff.png'),'Resolution',230);
savefig(tradeoff_fig,fullfile(result_folder,'m400_accuracy_tradeoff.fig'));

curve_fig = figure('Color','w','Position',[70 35 1350 880]);
curve_layout = tiledlayout(curve_fig,3,3, ...
    'TileSpacing','compact','Padding','compact');
colors = lines(numel(epsilon_values));
metric_titles = {'Control tracking error','Prediction error', ...
    'DAC average-tracking error'};
for row_i = 1:numel(round_values)
    for metric_i = 1:3
        ax = nexttile(curve_layout); hold(ax,'on');
        for col_i = 1:numel(epsilon_values)
            data = data_set{row_i,col_i};
            if metric_i == 1
                curve_t = data.t_set;
                curve_y = data.TrackingError_vector;
            elseif metric_i == 2
                curve_t = data.t_set;
                curve_y = data.prediction_error_norm_vector;
            else
                curve_t = data.t_set(1:numel(data.dac_tracking_error_set));
                curve_y = data.dac_tracking_error_set;
            end
            plot(ax,curve_t,curve_y,'Color',colors(col_i,:), ...
                'LineWidth',1.15,'DisplayName',sprintf( ...
                '\\epsilon=%.2f',epsilon_values(col_i)));
        end
        xlim(ax,[0 30]); grid(ax,'on'); box(ax,'on');
        title(ax,sprintf('R=%d: %s',round_values(row_i), ...
            metric_titles{metric_i}));
        if metric_i == 1, ylabel(ax,'Error norm'); end
        if row_i == numel(round_values), xlabel(ax,'Time (s)'); end
        if row_i == 1 && metric_i == 3
            legend(ax,'Location','best');
        end
    end
end
title(curve_layout,'M=400 IP-DAC accuracy across R and threshold');
exportgraphics(curve_fig,fullfile(result_folder, ...
    'm400_error_curves.png'),'Resolution',230);
savefig(curve_fig,fullfile(result_folder,'m400_error_curves.fig'));

disp(Summary);
end
