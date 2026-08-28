function Summary = run_single_poe_baselines(do_run, simulation_end_time)
%RUN_SINGLE_POE_BASELINES One controlled PoE comparison without Monte Carlo.
%
% Methods:
%   Exact, Offline-NBR (350 points/agent), TT-Local,
%   TT-NBR-PoE, TT-CEN-PoE, and converged TT-TP-PoE.
%
% All methods use the same formation task and initial-state seed.  The TT
% methods add one sample per agent every 0.1 s.  This script plots the same
% instantaneous overall tracking-error norm used in the shared code:
%       E(t) = norm(vartheta_all(t), 2).

if nargin < 1 || isempty(do_run)
    do_run = true;
end
if nargin < 2 || isempty(simulation_end_time)
    simulation_end_time = 30;
end

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));

result_folder = fullfile(repo_root, 'result', 'single_poe_baselines');
if ~exist(result_folder, 'dir')
    mkdir(result_folder);
end

old_end_time = getenv('CONTROL_SIM_END_TIME');
cleanup_end_time = onCleanup(@() setenv('CONTROL_SIM_END_TIME', old_end_time));
setenv('CONTROL_SIM_END_TIME', num2str(simulation_end_time, '%.15g'));

initial_seed = 42;
use_formation = true;
offline_points_per_agent = 350;

specs = [ ...
    struct('key','exact',       'label','Exact',           'runner','test', 'mode','exact',       'offline_points',0); ...
    struct('key','offline_nbr', 'label','Offline-NBR',     'runner','nbr',  'mode','poe_offline', 'offline_points',offline_points_per_agent); ...
    struct('key','tt_local',    'label','TT-Local',        'runner','test', 'mode','local',       'offline_points',0); ...
    struct('key','tt_nbr_poe',  'label','TT-NBR-PoE',      'runner','nbr',  'mode','poe',         'offline_points',0); ...
    struct('key','tt_cen_poe',  'label','TT-CEN-PoE',      'runner','cen',  'mode','poe',         'offline_points',0); ...
    struct('key','tt_tp_poe',   'label','TT-TP-PoE',       'runner','test', 'mode','poe_ac',      'offline_points',0)];

if do_run
    for method_idx = 1:numel(specs)
        spec = specs(method_idx);
        fprintf('\n[%d/%d] Running %s (seed=%d, T=%.1f s)\n', ...
            method_idx, numel(specs), spec.label, initial_seed, simulation_end_time);
        switch spec.runner
            case 'test'
                run_simulation_test_point(spec.mode, result_folder, spec.key, ...
                    use_formation, initial_seed, spec.offline_points);
            case 'nbr'
                run_simulation_nbr(spec.mode, result_folder, spec.key, ...
                    use_formation, initial_seed, spec.offline_points);
            case 'cen'
                run_simulation_cen(spec.mode, result_folder, spec.key, ...
                    use_formation, initial_seed, spec.offline_points);
            otherwise
                error('Unknown runner: %s', spec.runner);
        end
    end
end

method_count = numel(specs);
Method = strings(method_count,1);
InitialPointsPerAgent = zeros(method_count,1);
OnlineSamplesPerAgent = zeros(method_count,1);
FinalTrackingError = nan(method_count,1);
InitialTrackingError = nan(method_count,1);
FullTimeMeanTrackingError = nan(method_count,1);

fig = figure('Color','w','Position',[100 100 1050 620]);
ax = axes(fig);
hold(ax,'on');
colors = [0.0000 0.0000 0.0000; ...
          0.8500 0.3250 0.0980; ...
          0.0000 0.4470 0.7410; ...
          0.4940 0.1840 0.5560; ...
          0.4660 0.6740 0.1880; ...
          0.3010 0.7450 0.9330];
line_styles = {'-','--','-.',':','-','--'};

for method_idx = 1:method_count
    spec = specs(method_idx);
    result_path = fullfile(result_folder, [spec.key '.mat']);
    if ~isfile(result_path)
        error('Missing result file: %s', result_path);
    end
    data = load(result_path, 't_set', 'TrackingError_vector', ...
        'online_trigger_count');

    t = data.t_set(:);
    tracking = data.TrackingError_vector(:);
    if numel(t) ~= numel(tracking) || any(~isfinite(tracking))
        error('%s produced an invalid tracking-error curve.', spec.label);
    end

    semilogy(ax, t, max(tracking,eps), ...
        'LineStyle',line_styles{method_idx}, ...
        'LineWidth',1.8, ...
        'Color',colors(method_idx,:), 'DisplayName',spec.label);
    Method(method_idx) = spec.label;
    InitialPointsPerAgent(method_idx) = spec.offline_points;
    if isfield(data,'online_trigger_count') && ~isempty(data.online_trigger_count)
        OnlineSamplesPerAgent(method_idx) = mean(data.online_trigger_count);
    end
    InitialTrackingError(method_idx) = tracking(1);
    FinalTrackingError(method_idx) = tracking(end);
    FullTimeMeanTrackingError(method_idx) = trapz(t,tracking) / (t(end)-t(1));
end

xlabel(ax,'Time (s)');
ylabel(ax,'Overall tracking error ||\vartheta(t)||_2');
title(ax,'Single-run PoE baseline comparison');
legend(ax,'Location','best');
grid(ax,'on');
box(ax,'on');
set(ax,'YScale','log');
xlim(ax,[0 simulation_end_time]);
ylim(ax,[0.06 0.30]);
exportgraphics(fig, fullfile(result_folder,'tracking_error_single_run.png'), ...
    'Resolution',220);
savefig(fig, fullfile(result_folder,'tracking_error_single_run.fig'));

if max(InitialTrackingError)-min(InitialTrackingError) > 1e-12
    error('Methods did not start from the same tracking error.');
end
Summary = table(Method, InitialPointsPerAgent, OnlineSamplesPerAgent, ...
    InitialTrackingError, FinalTrackingError, FullTimeMeanTrackingError);
writetable(Summary, fullfile(result_folder,'single_run_summary.csv'));
save(fullfile(result_folder,'single_run_summary.mat'),'Summary','specs', ...
    'simulation_end_time','initial_seed','use_formation');

fprintf('\nSingle-run PoE summary\n');
disp(Summary);
fprintf('Results: %s\n', result_folder);
end
