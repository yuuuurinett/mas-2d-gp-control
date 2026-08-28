%% plot_M_ablation_final.m
% Final M-ablation plots for IP-DAC.
%
% Layout: 4 figures (one per dataset), each with 2x2 subplots:
%   Top-left:  SMSE    Top-right: RMSE
%   Bot-left:  MSLL    Bot-right: Train Time (ms/pt)
%
% One version saved per dataset:
%   _full : full M range
%   _zoom : zoomed in to the flat region (last ~1/3 of M range)
%
% Error bar style:
%   - Shaded band (mean ± std) across full M range
%   - Explicit error bars every few M points for clarity
%
% Output: Result/Figures/M_ablation_final/without_MOE/

clear; clc; close all;

ProjectRoot = fileparts(mfilename('fullpath'));
cd(ProjectRoot);
addpath(genpath(ProjectRoot));

%% ==================== Configuration ====================

datasets    = {'KIN40K', 'POL', 'PUMADYN32NM', 'SARCOS'};
train_ratio = 0.4;
tr_tag      = round(train_ratio * 100);
mc_list     = 1:10;
n_mc        = numel(mc_list);

M_default   = 100:100:2000;
M_puma      = 100:20:500;

% Zoom x-ranges (tightly focused on flat/convergence region)
zoom_default = [1000, 2000];   % flat/convergence region after cropping M > 2000
zoom_puma    = [360,  500];    % last portion of 100:20:500

%% ==================== Plot style ====================

LineWidth      = 2.2;
MarkerSize     = 4.5;
ShadeAlpha     = 0.45;    % stronger shaded band for mean ± std

AxisFontSize   = 12;
LabelFontSize  = 13;
TitleFontSize  = 14;
SgtFontSize    = 16;
LegFontSize    = 10;

%% ==================== Figure versions ====================

versions(1).name         = 'without_MOE';
versions(1).methods      = {'poe','gpoe','bcm','rbcm'};
versions(1).labels       = {'POE','GPOE','BCM','RBCM'};
versions(1).markers      = {'none','none','none','none'};

% MOE is intentionally omitted from the final M-ablation figures.

%% ==================== Metrics ====================
% {field_in_mat,  y_label,               transform}
% transform: 'direct' or 'sqrt' (sqrt of smse gives NRMSE)
metric_defs = {
    'smse',               'SMSE',                  'direct';
    'rmse',               'RMSE',                  'direct';
    'msll',               'MSLL',                  'direct';
    't_train_per_point',  'Train Time (ms/pt)',     'direct';
    't_test_per_point',   'Test Time (ms/pt)',      'direct';
};
n_metrics = size(metric_defs, 1);

%% ==================== Output folder ====================

OutRoot = fullfile(ProjectRoot, 'Result', 'Figures', 'M_ablation_final');
if ~exist(OutRoot, 'dir'), mkdir(OutRoot); end

%% ==================== Main loop ====================

for vi = 1:numel(versions)
    ver      = versions(vi);
    nA       = numel(ver.methods);
    colors   = lines(nA);
    VerDir   = fullfile(OutRoot, ver.name);
    if ~exist(VerDir, 'dir'), mkdir(VerDir); end

    fprintf('\n[%s] Generating figures...\n', ver.name);

    for di = 1:numel(datasets)
        dname = datasets{di};

        if strcmp(dname, 'PUMADYN32NM')
            M_list    = M_puma;
            zoom_xr   = zoom_puma;
            xtk       = 100:100:500;
            xtk_zoom  = 400:20:500;
        else
            M_list    = M_default;
            zoom_xr   = zoom_default;
            xtk       = 500:500:2000;
            xtk_zoom  = 1000:200:2000;
        end

        nM = numel(M_list);

        %% ---------- Load data ----------
        % data(method, M_idx, metric_idx, mc_idx)
        data = NaN(nA, nM, n_metrics, n_mc);

        for ai = 1:nA
            method = ver.methods{ai};
            for mi = 1:nM
                M = M_list(mi);
                for ci = 1:n_mc
                    fname = sprintf('%s_M%d_tr%d_mc%d.mat', method, M, tr_tag, mc_list(ci));
                    fpath = fullfile(ProjectRoot, 'Result', 'Dataset', dname, fname);
                    if ~isfile(fpath), continue; end
                    try
                        S = load(fpath);
                        for ki = 1:n_metrics
                            field = metric_defs{ki,1};
                            if isfield(S, field)
                                v = S.(field);
                                if isnumeric(v) && isscalar(v) && isfinite(v)
                                    if strcmp(metric_defs{ki,3}, 'sqrt')
                                        v = sqrt(v);
                                    end
                                    data(ai, mi, ki, ci) = v;
                                end
                            end
                        end
                    catch
                    end
                end
            end
        end

        %% ---------- Filter outliers in Train Time ----------
        % Occasionally a single MC run has anomalously long train time
        % (e.g. 40x median) due to system load spikes. Remove those points
        % so they don't distort the mean/std curves.
        for ki = 1:n_metrics
            if contains(metric_defs{ki,2}, 'Time')
                for ai = 1:nA
                    for mi = 1:nM
                        slice = squeeze(data(ai, mi, ki, :));
                        med_v = median(slice, 'omitnan');
                        if isfinite(med_v) && med_v > 0
                            slice(slice > 3 * med_v) = NaN;
                            data(ai, mi, ki, :) = slice;
                        end
                    end
                end
            end
        end

        %% ---------- Pre-compute mean/std ----------
        y_mean = mean(data, 4, 'omitnan');  % [nA x nM x n_metrics]
        y_std  = std(data,  0, 4, 'omitnan');

        %% ---------- Generate full + zoom ----------
        view_names  = {'full', 'zoom'};
        view_xranges = {[min(M_list), max(M_list)],  zoom_xr};
        view_xticks  = {xtk,                          xtk_zoom};

        for view_i = 1:2
            vname   = view_names{view_i};
            x_range = view_xranges{view_i};
            x_ticks = view_xticks{view_i};

            fig = figure('Color','w', 'Units','centimeters', 'Position',[2,2,36,22]);
            tlo = tiledlayout(2, 3, 'TileSpacing','compact', 'Padding','compact');

            leg_handles = gobjects(nA, 1);

            for ki = 1:n_metrics
                ax = nexttile(tlo);
                hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');

                ylabel_str = metric_defs{ki, 2};

                % x_mask for y-limit computation
                x_mask = M_list >= x_range(1) & M_list <= x_range(2);

                for ai = 1:nA
                    ym = squeeze(y_mean(ai, :, ki));   % [1 x nM]
                    ys = squeeze(y_std(ai,  :, ki));

                    if all(isnan(ym)), continue; end

                    c = colors(ai,:);
                    mk = ver.markers{ai};

                    %-- 1. Shaded mean±std band (full M range, skip for Train Time) --%
                    if true
                        xb  = M_list(:)';
                        yup = (ym + ys);
                        ylo = (ym - ys);
                        vld = isfinite(yup) & isfinite(ylo);
                        if any(vld)
                            patch(ax, [xb(vld), fliplr(xb(vld))], ...
                                [yup(vld), fliplr(ylo(vld))], c, ...
                                'FaceAlpha', ShadeAlpha, 'EdgeColor','none', ...
                                'HandleVisibility','off');
                        end
                    end

                    %-- 2. Mean curve with sparse markers --%
                    mk_idx = choose_marker_idx(M_list, x_range);
                    h = plot(ax, M_list, ym, '-', ...
                        'Color', c, 'LineWidth', LineWidth, ...
                        'Marker', mk, 'MarkerIndices', mk_idx, ...
                        'MarkerSize', MarkerSize, 'MarkerFaceColor', c, ...
                        'DisplayName', ver.labels{ai});
                    if ki == 1
                        leg_handles(ai) = h;
                    end

                end

                %-- Axis formatting --%
                xlim(ax, x_range);
                xticks(ax, x_ticks);

                % Tight y-limits based on visible range only
                all_y = squeeze(data(:, x_mask, ki, :));
                ymin_v = min(all_y(:), [], 'omitnan');
                ymax_v = max(all_y(:), [], 'omitnan');
                if isfinite(ymin_v) && isfinite(ymax_v) && ymax_v > ymin_v
                    rng_y = ymax_v - ymin_v;
                    ylim(ax, [ymin_v - 0.12*rng_y,  ymax_v + 0.18*rng_y]);
                end

                % Tick format
                if contains(ylabel_str, 'Time')
                    ytickformat(ax, '%.3f');
                elseif contains(ylabel_str, 'SMSE') || contains(ylabel_str, 'RMSE')
                    ytickformat(ax, '%.3f');
                end

                set(ax, 'FontName','Times New Roman', 'FontSize', AxisFontSize, 'LineWidth',1.1);
                xlabel(ax, 'Number of Inducing Points, M', ...
                    'FontName','Times New Roman', 'FontSize', LabelFontSize);
                ylabel(ax, ylabel_str, 'Interpreter','none', ...
                    'FontName','Times New Roman', 'FontSize', LabelFontSize);
                title(ax, ylabel_str, 'FontName','Times New Roman', ...
                    'FontSize', TitleFontSize, 'FontWeight','bold');
            end

            % Legend in the sixth tile
            ax_leg = nexttile(tlo);
            axis(ax_leg, 'off');
            valid_h = leg_handles(isgraphics(leg_handles));
            valid_l = ver.labels(isgraphics(leg_handles));
            legend(ax_leg, valid_h, valid_l, 'Location','northwest', 'Box','off', ...
                'FontName','Times New Roman', 'FontSize', LegFontSize);

            % Super title
            ds_label = strrep(dname, 'PUMADYN32NM', 'PUMADYN32NM');
            sgtitle(tlo, sprintf('%s: M-ablation under IP-DAC  (%s view)', ds_label, vname), ...
                'FontName','Times New Roman', 'FontSize', SgtFontSize, 'FontWeight','bold');

            % Save
            base = sprintf('M_ablation_%s_%s_%s', dname, ver.name, vname);
            exportgraphics(fig, fullfile(VerDir, [base '.png']), 'Resolution', 300);
            exportgraphics(fig, fullfile(VerDir, [base '.pdf']), 'ContentType','vector');
            savefig(fig,           fullfile(VerDir, [base '.fig']));
            fprintf('  Saved: %s\n', base);
            close(fig);
        end

        fprintf('[%s] %s done.\n', ver.name, dname);
    end
end

fprintf('\nAll figures saved to:\n%s\n', OutRoot);

%% ==================== Local functions ====================

function idx = choose_marker_idx(M_list, x_range)
% Return marker indices for plot(), only within x_range, spaced ~every 5 pts.
    in_range = M_list >= x_range(1) & M_list <= x_range(2);
    all_idx  = find(in_range);
    if isempty(all_idx)
        idx = [];
        return;
    end
    step = max(1, round(numel(all_idx) / 8));   % ~8 markers visible
    idx  = all_idx(1:step:end);
    if idx(end) ~= all_idx(end)
        idx = [idx, all_idx(end)];
    end
end
