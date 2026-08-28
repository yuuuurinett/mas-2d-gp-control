function output_file = generate_mc10_results_markdown( ...
    result_folder,shared_result_folder)
%GENERATE_MC10_RESULTS_MARKDOWN Export the corrected per-agent MC10 results.
%
% All reported tracking quantities are first computed within each trial and
% then summarized across the ten seeds. Communication is agent-level.

if nargin < 1 || isempty(result_folder)
    repo_root = fileparts(fileparts(mfilename('fullpath')));
    result_folder = fullfile(repo_root,'result', ...
        'mc10_all_aggregation_m400_T30_per_agent');
end
if nargin < 2 || isempty(shared_result_folder)
    shared_result_folder = result_folder;
    if contains(lower(result_folder),'m600')
        candidate = strrep(result_folder,'m600','m400');
        if isfolder(candidate), shared_result_folder = candidate; end
    end
end

aggregations = {'poe','gpoe','moe','bcm','rbcm'};
aggregation_labels = {'POE','GPOE','MOE','BCM','RBCM'};
mode_keys = {'IP_DAC','IP_AC','TP_DAC','TP_AC','CEN','NBR', ...
    'Offline','Exact'};
mode_labels = { ...
    'IP-DAC (1 update / 0.01 s)', ...
    'IP-AC (PETC, sigma=0.3, R=10)', ...
    'TP-DAC (R=10 / 0.1 s)', ...
    'TP-AC (R=10 / 0.1 s)', ...
    'CEN','NBR','Offline Local','Exact'};
seed_values = 42:51;
configuration_file = fullfile(result_folder,'seed_042','poe_ip_dac.mat');
configuration = load(configuration_file,'NumInducingPoints');
num_inducing_points = configuration.NumInducingPoints;

n_aggregation = numel(aggregations);
n_mode = numel(mode_keys);
n_seed = numel(seed_values);
mean_tracking = nan(n_aggregation,n_mode,n_seed);
final_tracking = nan(n_aggregation,n_mode,n_seed);
max_after_10 = nan(n_aggregation,n_mode,n_seed);
broadcasts = nan(n_aggregation,n_mode,n_seed);

for aggregation_i = 1:n_aggregation
    aggregation = aggregations{aggregation_i};
    for seed_i = 1:n_seed
        seed_folder = fullfile(result_folder,sprintf( ...
            'seed_%03d',seed_values(seed_i)));
        files = { ...
            [aggregation '_ip_dac.mat'], ...
            [aggregation '_ip_ac.mat'], ...
            [aggregation '_tp_dac_fixedR10.mat'], ...
            [aggregation '_tp_ac_fixedR10.mat'], ...
            [aggregation '_cen.mat'], ...
            [aggregation '_nbr.mat'], ...
            'local_offline.mat', ...
            'exact.mat'};

        for mode_i = 1:n_mode
            result_file = fullfile(seed_folder,files{mode_i});
            if mode_i >= 3 && ~strcmp(result_folder,shared_result_folder)
                result_file = fullfile(shared_result_folder, ...
                    sprintf('seed_%03d',seed_values(seed_i)),files{mode_i});
            end
            if ~isfile(result_file)
                error('Missing MC10 result: %s',result_file);
            end
            d = load(result_file);
            if isfield(d,'vartheta_all_set') && ...
                    ~isempty(d.vartheta_all_set)
                tracking = sqrt(sum(d.vartheta_all_set.^2,1));
            else
                tracking = d.TrackingError_vector(:).';
            end
            mean_tracking(aggregation_i,mode_i,seed_i) = ...
                mean(tracking,'omitnan');
            final_tracking(aggregation_i,mode_i,seed_i) = tracking(end);
            max_after_10(aggregation_i,mode_i,seed_i) = ...
                max(tracking(d.t_set >= 10),[],'omitnan');

            switch mode_i
                case 1
                    broadcasts(aggregation_i,mode_i,seed_i) = ...
                        d.dac_broadcasts_per_agent;
                case 2
                    broadcasts(aggregation_i,mode_i,seed_i) = ...
                        d.ac_broadcasts_per_agent;
                case {3,4}
                    broadcasts(aggregation_i,mode_i,seed_i) = ...
                        mean(d.tp_broadcast_count_per_agent,'omitnan');
                case 6
                    % NBR exchanges its current local prediction at every
                    % 0.01 s control interval: 3000 broadcasts in 30 s.
                    broadcasts(aggregation_i,mode_i,seed_i) = ...
                        numel(d.t_set)-1;
            end
        end
    end
end

mean_of = @(x) mean(x,3,'omitnan');
std_of = @(x) std(x,0,3,'omitnan');
mean_tracking_mean = mean_of(mean_tracking);
mean_tracking_std = std_of(mean_tracking);
final_tracking_mean = mean_of(final_tracking);
final_tracking_std = std_of(final_tracking);
max_after_10_mean = mean_of(max_after_10);
max_after_10_std = std_of(max_after_10);
broadcast_mean = mean_of(broadcasts);
broadcast_std = std_of(broadcasts);

output_file = fullfile(result_folder,'mc10_results_summary.md');
fid = fopen(output_file,'w');
if fid < 0, error('Cannot open output file: %s',output_file); end
cleanup_file = onCleanup(@() fclose(fid));

fprintf(fid,'# Corrected MC10 control comparison\n\n');
fprintf(fid,['Configuration: 10 seeds (42--51), 30 s simulation, ' ...
    'online data = 300 samples per agent, Offline Local = 350 samples ' ...
    'per agent, inducing points M=%d. IP-AC, TP-DAC and TP-AC use ' ...
    'R=10.\n\n'],num_inducing_points);
if ~strcmp(result_folder,shared_result_folder)
    fprintf(fid,['Only IP-DAC/IP-AC were rerun for this inducing-point ' ...
        'setting. TP, CEN, NBR, Offline Local and Exact are reused from ' ...
        'the completed M400 common-baseline batch because they do not ' ...
        'depend on M.\n\n']);
end
fprintf(fid,['The primary advisor metric is ' ...
    '`max_{t >= 10 s} ||vartheta_all(t)||_2`. Values are Monte Carlo ' ...
    'mean +/- sample standard deviation.\n\n']);

fprintf(fid,'## Tracking error\n\n');
fprintf(fid,'| Method | IP-DAC | IP-AC | TP-DAC | TP-AC | CEN | NBR | Local | Exact |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for aggregation_i = 1:n_aggregation
    fprintf(fid,'| %s |',aggregation_labels{aggregation_i});
    for mode_i = 1:n_mode
        fprintf(fid,' %s |',format_tracking_value( ...
            max_after_10_mean(aggregation_i,mode_i), ...
            max_after_10_std(aggregation_i,mode_i)));
    end
    fprintf(fid,'\n');
end

fprintf(fid,'\nLocal denotes the offline local-GP baseline.\n\n');

fprintf(fid,'## Communication statistics\n\n');
fprintf(fid,['| Method | IP-DAC comm. / agent | IP-AC comm. / agent | ' ...
    'TP-DAC comm. / agent | TP-AC comm. / agent | ' ...
    'NBR comm. / agent |\n']);
fprintf(fid,'|---|---:|---:|---:|---:|---:|\n');
for aggregation_i = 1:n_aggregation
    fprintf(fid,'| %s |',aggregation_labels{aggregation_i});
    for mode_i = 1:4
        fprintf(fid,' %s |',format_value( ...
            broadcast_mean(aggregation_i,mode_i), ...
            broadcast_std(aggregation_i,mode_i),1));
    end
    fprintf(fid,' %s |\n',format_value( ...
        broadcast_mean(aggregation_i,6), ...
        broadcast_std(aggregation_i,6),1));
end

fprintf(fid,['\nAll communication values count agent-level packaged ' ...
    'broadcast instances; they are not per-inducing-point counts.\n\n']);

fprintf(fid,'## Communication-count definition\n\n');
fprintf(fid,['- IP-DAC and IP-AC: counts stored by the original ' ...
    'agent-level event-trigger implementation.\n']);
fprintf(fid,['- TP-DAC and TP-AC: 300 online snapshots x 10 fixed ' ...
    'consensus rounds = 3000 broadcasts per agent.\n']);
fprintf(fid,['- NBR: one scheduled neighbor broadcast per 0.01 s ' ...
    'control interval = 3000 broadcasts per agent.\n']);
fprintf(fid,['- CEN, Local and Exact: not reported as distributed ' ...
    'agent-level broadcast counts.\n']);

fprintf('Saved MC10 Markdown summary: %s\n',output_file);
end

function text = format_value(mean_value,std_value,decimal_places)
if ~isfinite(mean_value)
    text = '--';
    return;
end
format = sprintf('%%.%df +/- %%.%df',decimal_places,decimal_places);
text = sprintf(format,mean_value,std_value);
end

function text = format_tracking_value(mean_value,std_value)
if ~isfinite(mean_value)
    text = '--';
elseif std_value > 0 && std_value < 1e-5
    text = sprintf('%.5f +/- %.2e',mean_value,std_value);
else
    text = sprintf('%.5f +/- %.5f',mean_value,std_value);
end
end

% The code below is intentionally unreachable and retained only in older
% generated files; the current report uses the two compact tables above.
%{
fprintf(fid,'| Aggregation |');
for mode_i = 1:n_mode
    fprintf(fid,' %s |',strrep(mode_keys{mode_i},'_','-'));
end
fprintf(fid,'\n|---|');
for mode_i = 1:n_mode, fprintf(fid,'---:|'); end
fprintf(fid,'\n');
for aggregation_i = 1:n_aggregation
    fprintf(fid,'| %s |',aggregation_labels{aggregation_i});
    for mode_i = 1:n_mode
        fprintf(fid,' %s |',format_tracking_value( ...
            max_after_10_mean(aggregation_i,mode_i), ...
            max_after_10_std(aggregation_i,mode_i)));
    end
    fprintf(fid,'\n');
end

fprintf(fid,'\n## Detailed metrics and communication\n\n');
fprintf(fid,['`MeanTrackingError` is the full-trajectory time mean in ' ...
    'each trial. `FinalTrackingError` is the value at 30 s. ' ...
    '`BroadcastsPerAgent` counts agent-level broadcast instances.\n\n']);

for aggregation_i = 1:n_aggregation
    fprintf(fid,'### %s\n\n',aggregation_labels{aggregation_i});
    fprintf(fid,['| Method | MeanTrackingError | FinalTrackingError | ' ...
        'MaxTrackingErrorAfter10s | BroadcastsPerAgent |\n']);
    fprintf(fid,'|---|---:|---:|---:|---:|\n');
    for mode_i = 1:n_mode
        broadcast_text = format_value(broadcast_mean(aggregation_i,mode_i), ...
            broadcast_std(aggregation_i,mode_i),1);
        fprintf(fid,'| %s | %s | %s | %s | %s |\n', ...
            mode_labels{mode_i}, ...
            format_tracking_value( ...
                mean_tracking_mean(aggregation_i,mode_i), ...
                mean_tracking_std(aggregation_i,mode_i)), ...
            format_tracking_value( ...
                final_tracking_mean(aggregation_i,mode_i), ...
                final_tracking_std(aggregation_i,mode_i)), ...
            format_tracking_value( ...
                max_after_10_mean(aggregation_i,mode_i), ...
                max_after_10_std(aggregation_i,mode_i)),broadcast_text);
    end
    fprintf(fid,'\n');
end

fprintf(fid,'## Communication-count definition\n\n');
fprintf(fid,['- IP-DAC and IP-AC: counts stored by the original ' ...
    'agent-level event-trigger implementation.\n']);
fprintf(fid,['- TP-DAC and TP-AC: 300 online snapshots x 10 fixed ' ...
    'consensus rounds = 3000 broadcasts per agent.\n']);
fprintf(fid,['- NBR: one scheduled neighbor broadcast per 0.01 s ' ...
    'control interval = 3000 broadcasts per agent.\n']);
fprintf(fid,['- CEN, Offline Local and Exact: not reported as distributed ' ...
    'agent-level broadcast counts.\n']);

fprintf('Saved MC10 Markdown summary: %s\n',output_file);
end

function text = format_value(mean_value,std_value,decimal_places)
if ~isfinite(mean_value)
    text = '--';
    return;
end
format = sprintf('%%.%df +/- %%.%df',decimal_places,decimal_places);
text = sprintf(format,mean_value,std_value);
end

function text = format_tracking_value(mean_value,std_value)
if ~isfinite(mean_value)
    text = '--';
elseif std_value > 0 && std_value < 1e-5
    text = sprintf('%.5f +/- %.2e',mean_value,std_value);
else
    text = sprintf('%.5f +/- %.5f',mean_value,std_value);
end
%}
