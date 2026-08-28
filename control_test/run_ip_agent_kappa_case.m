function run_ip_agent_kappa_case(kappa_value)
%RUN_IP_AGENT_KAPPA_CASE One M=400, R=1, epsilon=0.1 Kappa_P case.

validateattributes(kappa_value,{'numeric'}, ...
    {'scalar','positive','finite'});
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result','ip_agent_kappa_sweep');
if ~isfolder(result_folder), mkdir(result_folder); end

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_DAC_STEPS_PER_PHYSICAL_STEP','CONTROL_IP_PROJECTION_DIAGNOSTICS', ...
    'CONTROL_IP_INDUCING_POINT_FILE','CONTROL_DAC_TRIGGER_EPSILON', ...
    'CONTROL_DAC_KAPPA_P'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment(env_names,old_values)); %#ok<NASGU>
setenv(env_names{1},'30');
setenv(env_names{2},'400');
setenv(env_names{3},'1');
setenv(env_names{4},'1');
setenv(env_names{5},'');
setenv(env_names{6},'0.1');
setenv(env_names{7},num2str(kappa_value,'%.15g'));

kappa_tag = strrep(sprintf('%.3g',kappa_value),'.','p');
result_name = ['poe_M400_R1_eps_0p1_kappa_',kappa_tag,'_T30'];
run_simulation_inducing_point('poe',result_folder,result_name, ...
    false,[],[],[],42);
end

function restore_environment(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
