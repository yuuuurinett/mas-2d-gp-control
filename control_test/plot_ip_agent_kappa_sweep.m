function Summary = plot_ip_agent_kappa_sweep()
%PLOT_IP_AGENT_KAPPA_SWEEP Diagnose DAC speed versus IP approximation.

repo_root = fileparts(fileparts(mfilename('fullpath')));
result_folder = fullfile(repo_root,'result','ip_agent_kappa_sweep');
kappa_values = [0.5 1 2 5];
n = numel(kappa_values);
KappaP = kappa_values(:);
MeanBroadcastRate = nan(n,1); Last5sBroadcastRate = nan(n,1);
MeanTrackingError = nan(n,1); MeanIPDACError = nan(n,1);
MeanIPCENError = nan(n,1); MeanDirectError = nan(n,1);
MeanDACTrackingError = nan(n,1); LateDACTrackingError = nan(n,1);
data_set = cell(n,1); rates = cell(n,1);

for idx = 1:n
    tag = strrep(sprintf('%.3g',kappa_values(idx)),'.','p');
    data = load(fullfile(result_folder, ...
        ['poe_M400_R1_eps_0p1_kappa_',tag,'_T30.mat']));
    data_set{idx} = data;
    events = double(data.dac_inner_trigger_instance_set);
    event_t = data.t_set(1:size(events,3));
    rate = nan(1,30);
    for second_i = 1:30
        mask = event_t >= second_i-1 & event_t < second_i;
        rate(second_i) = sum(events(:,:,mask),'all')/size(events,1);
    end
    rates{idx} = rate;
    MeanBroadcastRate(idx) = mean(rate);
    Last5sBroadcastRate(idx) = mean(rate(end-4:end));
    MeanTrackingError(idx) = trapz( ...
        data.t_set,data.TrackingError_vector)/data.t_set(end);
    valid = isfinite(data.direct_prediction_error_norm_vector) & ...
        isfinite(data.ip_cen_prediction_error_norm_vector) & ...
        isfinite(data.prediction_error_norm_vector);
    MeanDirectError(idx) = mean( ...
        data.direct_prediction_error_norm_vector(valid));
    MeanIPCENError(idx) = mean( ...
        data.ip_cen_prediction_error_norm_vector(valid));
    MeanIPDACError(idx) = mean( ...
        data.prediction_error_norm_vector(valid));
    MeanDACTrackingError(idx) = mean( ...
        data.dac_tracking_error_set,'omitnan');
    LateDACTrackingError(idx) = mean( ...
        data.dac_tracking_error_set(event_t>=25),'omitnan');
end

Summary = table(KappaP,MeanBroadcastRate,Last5sBroadcastRate, ...
    MeanTrackingError,MeanDirectError,MeanIPCENError,MeanIPDACError, ...
    MeanDACTrackingError,LateDACTrackingError);
writetable(Summary,fullfile(result_folder,'summary.csv'));
save(fullfile(result_folder,'summary.mat'),'Summary');

fig = figure('Color','w','Position',[80 70 1250 760]);
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
colors = lines(n);
ax = nexttile(layout); hold(ax,'on');
for idx = 1:n
    plot(ax,(0:29)+0.5,rates{idx},'Color',colors(idx,:), ...
        'LineWidth',1.4,'DisplayName',sprintf('\\kappa_P=%.1f',KappaP(idx)));
end
xlabel(ax,'Time (s)'); ylabel(ax,'Broadcasts/agent/s');
title(ax,'Packaged communication'); grid(ax,'on'); box(ax,'on'); legend(ax);

ax = nexttile(layout); hold(ax,'on');
plot(ax,KappaP,MeanDACTrackingError,'o-','LineWidth',1.7, ...
    'DisplayName','Full trajectory');
plot(ax,KappaP,LateDACTrackingError,'s-','LineWidth',1.7, ...
    'DisplayName','Last 5 s');
xlabel(ax,'\\kappa_P'); ylabel(ax,'DAC average-tracking error');
title(ax,'Tracking current mean P'); grid(ax,'on'); box(ax,'on'); legend(ax);

ax = nexttile(layout); hold(ax,'on');
plot(ax,KappaP,MeanDirectError,'o-','LineWidth',1.7,'DisplayName','Direct PoE');
plot(ax,KappaP,MeanIPCENError,'s-','LineWidth',1.7,'DisplayName','IP-CEN');
plot(ax,KappaP,MeanIPDACError,'d-','LineWidth',1.7,'DisplayName','IP-DAC');
xlabel(ax,'\\kappa_P'); ylabel(ax,'Mean prediction error');
title(ax,'Prediction bottleneck'); grid(ax,'on'); box(ax,'on'); legend(ax);

ax = nexttile(layout); hold(ax,'on');
yyaxis(ax,'left');
plot(ax,KappaP,MeanBroadcastRate,'o-','LineWidth',1.7);
ylabel(ax,'Broadcasts/agent/s');
yyaxis(ax,'right');
plot(ax,KappaP,MeanTrackingError,'s-','LineWidth',1.7);
ylabel(ax,'Mean control tracking error');
xlabel(ax,'\\kappa_P'); title(ax,'Communication-control tradeoff');
grid(ax,'on'); box(ax,'on');
title(layout,'M=400, R=1, \\epsilon=0.1: DAC gain diagnosis');
exportgraphics(fig,fullfile(result_folder,'kappa_diagnosis.png'),'Resolution',230);
savefig(fig,fullfile(result_folder,'kappa_diagnosis.fig'));

best_idx = find(MeanDACTrackingError == min(MeanDACTrackingError),1);
data = data_set{best_idx};
valid = isfinite(data.direct_prediction_error_norm_vector) & ...
    isfinite(data.ip_cen_prediction_error_norm_vector) & ...
    isfinite(data.prediction_error_norm_vector);
comparison_fig = figure('Color','w','Position',[100 100 1150 500]);
ax = axes(comparison_fig); hold(ax,'on');
plot(ax,data.t_set(valid),data.direct_prediction_error_norm_vector(valid), ...
    'LineWidth',1.5,'DisplayName','Direct PoE');
plot(ax,data.t_set(valid),data.ip_cen_prediction_error_norm_vector(valid), ...
    'LineWidth',1.5,'DisplayName','IP-CEN');
plot(ax,data.t_set(valid),data.prediction_error_norm_vector(valid), ...
    'LineWidth',1.5,'DisplayName','IP-DAC');
xlabel(ax,'Time (s)'); ylabel(ax,'Prediction error norm');
title(ax,sprintf('Prediction comparison at \\kappa_P=%.1f',KappaP(best_idx)));
grid(ax,'on'); box(ax,'on'); legend(ax,'Location','best');
exportgraphics(comparison_fig,fullfile(result_folder, ...
    'best_kappa_prediction_comparison.png'),'Resolution',230);
savefig(comparison_fig,fullfile(result_folder, ...
    'best_kappa_prediction_comparison.fig'));
disp(Summary);
end
