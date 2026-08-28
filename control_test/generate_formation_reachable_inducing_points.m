function [output_file,position_domain_scale] = ...
    generate_formation_reachable_inducing_points(output_file,margin_fraction)
%GENERATE_FORMATION_REACHABLE_INDUCING_POINTS Known-reference IP domain.
% The domain uses only the prescribed leader/formation trajectories, not
% simulated states or unknown-dynamics observations.  M and all GP
% hyperparameters remain unchanged: 100 position pairs x 6 velocity pairs.

repo_root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || isempty(output_file)
    output_file = fullfile(repo_root,'result', ...
        'inducing_points_position100_velocity6_formation_reachable.mat');
end
if nargin < 2 || isempty(margin_fraction), margin_fraction = 0.10; end

reference_time = 0:0.01:20;
[leader_state,~,~] = Manipulator_2D_2DoF_LeaderDynamics(reference_time,1);
agent_quantity = 6;
desired_positions = nan(2,agent_quantity,numel(reference_time));
for agent_i = 1:agent_quantity
    formation_state = Manipulator_2D_2DoF_RelativePositionDynamics( ...
        reference_time,agent_i,agent_quantity);
    desired_position_i = leader_state(1:2,:) + formation_state(1:2,:);
    desired_positions(:,agent_i,:) = reshape(desired_position_i, ...
        2,1,numel(reference_time));
end

maximum_reference_position = max(abs(desired_positions),[],'all');
position_domain_scale = (1+margin_fraction)*maximum_reference_position;
generate_reference_velocity_inducing_points(output_file,100,6, ...
    position_domain_scale);

sampling_source = 'known leader plus formation reference, no rollout data';
save(output_file,'margin_fraction','maximum_reference_position', ...
    'sampling_source','-append');
fprintf('Formation-reference position domain: [%.3f, %.3f].\n', ...
    -position_domain_scale,position_domain_scale);
end
