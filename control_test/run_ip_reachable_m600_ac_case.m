function run_ip_reachable_m600_ac_case( ...
    simulation_end_time,consensus_tol,consensus_step,trigger_sigma, ...
    stop_mode)
%RUN_IP_REACHABLE_M600_AC_CASE Packaged-trigger static AC comparison.

if nargin < 1 || isempty(simulation_end_time), simulation_end_time = 40; end
if nargin < 2 || isempty(consensus_tol), consensus_tol = 0.05; end
if nargin < 3 || isempty(consensus_step), consensus_step = 0.1; end
if nargin < 4 || isempty(trigger_sigma), trigger_sigma = 0.5; end
if nargin < 5 || isempty(stop_mode), stop_mode = 'tracking'; end
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result','ip_reachable_m600');
point_file = fullfile(result_folder,'reachable_uniform_M600.mat');

env_names = {'CONTROL_SIM_END_TIME','CONTROL_IP_NUM_INDUCING_POINTS', ...
    'CONTROL_IP_PROJECTION_DIAGNOSTICS','CONTROL_IP_INDUCING_POINT_FILE', ...
    'CONTROL_DAC_KAPPA_P','CONTROL_AC_TRIGGER_SIGMA', ...
    'CONTROL_CONSENSUS_TOL','CONTROL_AC_CONSENSUS_STEP'};
env_names{end+1} = 'CONTROL_AC_STOP_MODE';
old_values = cellfun(@getenv,env_names,'UniformOutput',false);
cleanup_env = onCleanup(@() restore_environment(env_names,old_values)); %#ok<NASGU>
setenv(env_names{1},num2str(simulation_end_time,'%.15g'));
setenv(env_names{2},'600'); setenv(env_names{3},'0');
setenv(env_names{4},point_file); setenv(env_names{5},'1');
setenv(env_names{6},num2str(trigger_sigma,'%.15g'));
setenv(env_names{7},num2str(consensus_tol,'%.15g'));
setenv(env_names{8},num2str(consensus_step,'%.15g'));
setenv(env_names{9},stop_mode);
tol_tag = strrep(sprintf('%.3g',consensus_tol),'.','p');
step_tag = strrep(sprintf('%.3g',consensus_step),'.','p');
sigma_tag = strrep(sprintf('%.3g',trigger_sigma),'.','p');
run_simulation_inducing_point('poe_ac',result_folder, ...
    sprintf(['poe_ac_M600_reachable_packaged_sigma_%s_', ...
    'tol_%s_step_%s_stop_%s_T%g'],sigma_tag,tol_tag,step_tag, ...
    stop_mode,simulation_end_time), ...
    false,[],[],[],42);
end

function restore_environment(names,values)
for idx = 1:numel(names), setenv(names{idx},values{idx}); end
end
