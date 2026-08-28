function [TrackingTable,CommunicationTable] = ...
    summarize_current_seed(seed,use_formation)
%SUMMARIZE_CURRENT_SEED Summarize one completed formal comparison seed.
if nargin < 1, seed = 1; end
if nargin < 2 || isempty(use_formation), use_formation = false; end

repo_root = fileparts(fileparts(mfilename('fullpath')));
result_name = 'mc50_all_modes_T20_U0p4_D0p1_M600_structured_R10';
if use_formation, result_name = [result_name '_formation']; end
result_root = fullfile(repo_root,'result',result_name, ...
    sprintf('seed_%03d',seed));
methods = {'poe','gpoe','moe','bcm','rbcm'};
mode_labels = {'IP_DAC','IP_AC','TP_DAC','TP_AC','CEN','NBR', ...
    'Offline','Exact'};
suffixes = {'ip_dac','ip_ac','tp_dac','tp_ac','cen','nbr'};

tracking = nan(numel(methods),numel(mode_labels));
communication = nan(numel(methods),numel(mode_labels));
offline = load(fullfile(result_root,'offline.mat'));
exact = load(fullfile(result_root,'exact.mat'));
offline_metric = tracking_metric(offline);
exact_metric = tracking_metric(exact);

for method_idx = 1:numel(methods)
    method = methods{method_idx};
    for mode_idx = 1:numel(suffixes)
        d = load(fullfile(result_root, ...
            [method '_' suffixes{mode_idx} '.mat']));
        tracking(method_idx,mode_idx) = tracking_metric(d);
        switch mode_idx
            case 1
                communication(method_idx,mode_idx) = ...
                    mean(d.dac_broadcast_event_count_per_agent);
            case 2
                communication(method_idx,mode_idx) = ...
                    mean(d.ac_event_count_per_agent);
            case {3,4}
                communication(method_idx,mode_idx) = ...
                    mean(d.tp_broadcast_count_per_agent);
            case 6
                % NBR evaluates neighbor GPs at every 0.01 s physical step.
                communication(method_idx,mode_idx) = numel(d.t_set)-1;
        end
    end
    tracking(method_idx,7:8) = [offline_metric,exact_metric];
end

row_names = upper(string(methods));
TrackingTable = array2table(tracking,'VariableNames',mode_labels, ...
    'RowNames',cellstr(row_names));
CommunicationTable = array2table(communication, ...
    'VariableNames',mode_labels,'RowNames',cellstr(row_names));

disp('Tracking metric: max_{t >= 10 s} ||vartheta(t)||_2');
disp(TrackingTable);
disp('Agent-level broadcasts per agent over the 20 s trajectory');
disp(CommunicationTable);
end

function value = tracking_metric(d)
if isfield(d,'vartheta_all_set') && ~isempty(d.vartheta_all_set)
    curve = sqrt(sum(d.vartheta_all_set.^2,1));
else
    curve = d.TrackingError_vector;
end
mask = d.t_set >= 10;
value = max(curve(mask),[],'omitnan');
end
