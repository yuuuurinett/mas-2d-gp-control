function output_file = generate_reference_velocity_inducing_points( ...
    output_file,position_count,velocity_count,position_domain_scale)
%GENERATE_REFERENCE_VELOCITY_INDUCING_POINTS Structured position/velocity IP set.
% Position pairs uniformly cover [-1.5,1.5]^2. Velocity pairs are k-means
% representatives of the known leader reference trajectory, not samples
% from a control-test realization.

repo_root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || isempty(output_file)
    output_file = fullfile(repo_root,'result', ...
        'inducing_points_position100_velocity6_reference_seed42.mat');
end
if nargin < 2 || isempty(position_count), position_count = 100; end
if nargin < 3 || isempty(velocity_count), velocity_count = 6; end
if nargin < 4 || isempty(position_domain_scale), position_domain_scale = 1.5; end

reference_time = 0:0.01:20;
[leader_state,~,~] = Manipulator_2D_2DoF_LeaderDynamics( ...
    reference_time,1);

rng_before = rng;
cleanup_rng = onCleanup(@() rng(rng_before));
rng(42,'twister');
[~,velocity_centers] = kmeans(leader_state(3:4,:).',velocity_count, ...
    'Replicates',20,'MaxIter',1000);
velocity_centers = sortrows(velocity_centers,1);

generate_position_velocity_inducing_points(output_file,position_count, ...
    velocity_count,position_domain_scale,42, ...
    velocity_centers.');
reference_source = 'known-leader-reference-trajectory';
save(output_file,'reference_time','velocity_centers', ...
    'reference_source','position_domain_scale','-append');
end
