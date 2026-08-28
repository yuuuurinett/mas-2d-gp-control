function Stats = run_ip_dac_trigger_30s(do_run)
%RUN_IP_DAC_TRIGGER_30S Long-horizon trigger diagnostic with random IPs.

if nargin < 1 || isempty(do_run), do_run = true; end
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result','ip_dac_trigger_30s');
if ~isfolder(result_folder), mkdir(result_folder); end
result_name = 'poe_random_M600_R10_T30';
result_file = fullfile(result_folder,[result_name,'.mat']);

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_DAC_STEPS_PER_PHYSICAL_STEP','CONTROL_IP_PROJECTION_DIAGNOSTICS', ...
    'CONTROL_IP_INDUCING_POINT_FILE'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment(env_names,old_values)); %#ok<NASGU>
setenv(env_names{1},'30');
setenv(env_names{2},'600');
setenv(env_names{3},'10');
setenv(env_names{4},'0');
setenv(env_names{5},''); % Explicitly restore uniform random domain sampling.

if do_run
    run_simulation_inducing_point('poe',result_folder,result_name, ...
        false,[],[],[],42);
end
Stats = plot_ip_dac_interval_convergence(result_file);
end

function restore_environment(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
