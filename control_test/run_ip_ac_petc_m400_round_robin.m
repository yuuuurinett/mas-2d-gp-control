function result_file = run_ip_ac_petc_m400_round_robin(do_run)
%RUN_IP_AC_PETC_M400_ROUND_ROBIN Legacy entry point for the IP-AC detector.
% M=400, sigma=0.3, consensus gain=0.2, ten communication rounds per
% 0.1 s GP refresh, with one local sample added by every agent per refresh.

if nargin < 1 || isempty(do_run), do_run = true; end

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result', ...
    'ip_ac_petc_m400_all_agents_T30');
if ~isfolder(result_folder), mkdir(result_folder); end
result_name = 'poe_ac_petc_sigma03_gain02_R10_M400_online300_per_agent_T30';
result_file = fullfile(result_folder,[result_name '.mat']);

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_IP_INDUCING_POINT_FILE','CONTROL_AC_ITERATION_POLICY', ...
    'CONTROL_AC_FIXED_ITERATIONS','CONTROL_AC_BROADCAST_TRIGGER', ...
    'CONTROL_AC_PERIODIC_SIGMA','CONTROL_AC_TRIGGER_DIAGNOSTICS', ...
    'CONTROL_ONLINE_AGENT_POLICY','CONTROL_AC_CONSENSUS_STEP'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment_local( ...
    env_names,old_values)); %#ok<NASGU>

setenv(env_names{1},'30');
setenv(env_names{2},'400');
setenv(env_names{3},'');             % random inducing locations, fixed by seed
setenv(env_names{4},'fixed');
setenv(env_names{5},'10');
setenv(env_names{6},'petc_cache');
setenv(env_names{7},'0.3');
setenv(env_names{8},'1');            % save both sides of the detector
setenv(env_names{9},'all_agents');   % 300 samples in each local GP over 30 s
setenv(env_names{10},'0.2');         % communication-saving R10 compromise

if do_run || ~isfile(result_file)
    run_simulation_inducing_point('poe_ac',result_folder,result_name, ...
        false,[],[],[],42);
end

plot_ip_ac_advisor_detector(result_file,4,10.0);
end

function restore_environment_local(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
