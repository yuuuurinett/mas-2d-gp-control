function run_ip_agent_m400_grid_case(round_count,epsilon_value)
%RUN_IP_AGENT_M400_GRID_CASE One 30-s M=400 packaged-trigger grid case.

% DAC inner rounds subdivide one physical step, so the simulator uses
% t_step/round_count for each DAC update.

validateattributes(round_count,{'numeric'}, ...
    {'scalar','integer','positive','finite'});
validateattributes(epsilon_value,{'numeric'}, ...
    {'scalar','positive','finite'});

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result','ip_agent_m400_grid');
if ~isfolder(result_folder), mkdir(result_folder); end

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_DAC_STEPS_PER_PHYSICAL_STEP','CONTROL_IP_PROJECTION_DIAGNOSTICS', ...
    'CONTROL_IP_INDUCING_POINT_FILE','CONTROL_DAC_TRIGGER_EPSILON'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment(env_names,old_values)); %#ok<NASGU>
setenv(env_names{1},'30');
setenv(env_names{2},'400');
setenv(env_names{3},num2str(round_count));
setenv(env_names{4},'0');
setenv(env_names{5},'');
setenv(env_names{6},num2str(epsilon_value,'%.15g'));

epsilon_tag = strrep(sprintf('%.3g',epsilon_value),'.','p');
result_name = sprintf('poe_M400_R%d_eps_%s_T30',round_count,epsilon_tag);
run_simulation_inducing_point('poe',result_folder,result_name, ...
    false,[],[],[],42);
end

function restore_environment(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
