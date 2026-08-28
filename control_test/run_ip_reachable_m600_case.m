function run_ip_reachable_m600_case(epsilon_value,simulation_end_time)
%RUN_IP_REACHABLE_M600_CASE M=600 uniform reachable-box IP-DAC case.

validateattributes(epsilon_value,{'numeric'}, ...
    {'scalar','positive','finite'});
if nargin < 2 || isempty(simulation_end_time), simulation_end_time = 30; end
validateattributes(simulation_end_time,{'numeric'}, ...
    {'scalar','positive','finite'});
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result','ip_reachable_m600');
if ~isfolder(result_folder), mkdir(result_folder); end

point_source = load(fullfile(repo_root,'result','ip_reachable_uniform', ...
    'reachable_uniform_M400.mat'),'sampling_min','sampling_max');
rng(42,'twister');
M = 600;
InducingPoints_Coordinates = point_source.sampling_min + ...
    (point_source.sampling_max-point_source.sampling_min).*rand(4,M);
point_file = fullfile(result_folder,'reachable_uniform_M600.mat');
sampling_min = point_source.sampling_min; %#ok<NASGU>
sampling_max = point_source.sampling_max; %#ok<NASGU>
save(point_file,'InducingPoints_Coordinates','sampling_min','sampling_max');

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_DAC_STEPS_PER_PHYSICAL_STEP','CONTROL_IP_PROJECTION_DIAGNOSTICS', ...
    'CONTROL_IP_INDUCING_POINT_FILE','CONTROL_DAC_TRIGGER_EPSILON', ...
    'CONTROL_DAC_KAPPA_P'};
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment(env_names,old_values)); %#ok<NASGU>
setenv(env_names{1},num2str(simulation_end_time,'%.15g'));
setenv(env_names{2},'600');
setenv(env_names{3},'1'); setenv(env_names{4},'1');
setenv(env_names{5},point_file);
setenv(env_names{6},num2str(epsilon_value,'%.15g'));
setenv(env_names{7},'1');

epsilon_tag = strrep(sprintf('%.3g',epsilon_value),'.','p');
run_simulation_inducing_point('poe',result_folder, ...
    sprintf('poe_M600_reachable_R1_eps_%s_kappa_1_T%g', ...
    epsilon_tag,simulation_end_time), ...
    false,[],[],[],42);
end

function restore_environment(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
