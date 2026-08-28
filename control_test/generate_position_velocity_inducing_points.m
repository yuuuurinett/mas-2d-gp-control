function output_file = generate_position_velocity_inducing_points( ...
    output_file,position_count,velocity_count,domain_scale,seed, ...
    velocity_points_override)
%GENERATE_POSITION_VELOCITY_INDUCING_POINTS Structured 4-D inducing set.
% Use 100 well-spread two-dimensional positions and six representative
% two-dimensional velocities, then take their Cartesian product.

if nargin < 1 || isempty(output_file)
    repo_root = fileparts(fileparts(mfilename('fullpath')));
    output_file = fullfile(repo_root,'result', ...
        'inducing_points_position100_velocity6_seed42.mat');
end
if nargin < 2 || isempty(position_count), position_count = 100; end
if nargin < 3 || isempty(velocity_count), velocity_count = 6; end
if nargin < 4 || isempty(domain_scale), domain_scale = 1.5; end
if nargin < 5 || isempty(seed), seed = 42; end
if nargin < 6, velocity_points_override = []; end

if position_count*velocity_count ~= 600
    warning('The requested structured set contains %d points, not 600.', ...
        position_count*velocity_count);
end

rng_state_before_sampling = rng;
cleanup_rng = onCleanup(@() rng(rng_state_before_sampling)); %#ok<NASGU>
rng(seed,'twister');

grid_side = round(sqrt(position_count));
if grid_side^2 == position_count
    position_axis = linspace(-domain_scale,domain_scale,grid_side);
    [position_q1,position_q2] = meshgrid(position_axis,position_axis);
    position_points = [position_q1(:)';position_q2(:)'];
else
    % A non-square position count (for example 120 or 150) cannot form a
    % square tensor grid.  Use a reproducible maximin Latin hypercube so
    % both position coordinates remain equally and uniformly covered.
    % Plain random points would confound the allocation comparison with
    % accidental clustering and holes in the position domain.
    position_unit = lhsdesign(position_count,2, ...
        'Criterion','maximin','Iterations',100);
    position_points = (2*domain_scale*position_unit-domain_scale).';
end

if isempty(velocity_points_override)
    % Use the same domain as the position coordinates.  A Latin-hypercube
    % construction gives one representative in every velocity stratum;
    % only the number of velocity samples is reduced to six.
    velocity_fraction = ((1:velocity_count)-0.5)/velocity_count;
    velocity_permutation = randperm(velocity_count);
    velocity_points = domain_scale*(2*[ ...
        velocity_fraction; ...
        velocity_fraction(velocity_permutation)]-1);
else
    velocity_points = velocity_points_override;
    if ~isequal(size(velocity_points),[2,velocity_count])
        error('velocity_points_override must be 2-by-%d.',velocity_count);
    end
end
position_indices = repelem(1:position_count,velocity_count);
velocity_indices = repmat(1:velocity_count,1,position_count);
InducingPoints_Coordinates = [ ...
    position_points(:,position_indices); ...
    velocity_points(:,velocity_indices)];

sampling_scheme = 'position-grid-times-representative-velocities';
sampling_seed = seed;
output_folder = fileparts(output_file);
if ~isfolder(output_folder), mkdir(output_folder); end
save(output_file,'InducingPoints_Coordinates','position_points', ...
    'velocity_points','position_count','velocity_count','domain_scale', ...
    'sampling_scheme','sampling_seed');

fig = figure('Visible','off','Color','w','Position',[100 100 900 380]);
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(layout); scatter(ax1,position_points(1,:), ...
    position_points(2,:),18,'filled');
axis(ax1,'equal'); grid(ax1,'on'); box(ax1,'on');
xlim(ax1,domain_scale*[-1 1]); ylim(ax1,domain_scale*[-1 1]);
xlabel(ax1,'q_1'); ylabel(ax1,'q_2');
title(ax1,sprintf('Position samples (%d)',position_count), ...
    'FontWeight','normal');
ax2 = nexttile(layout); scatter(ax2,velocity_points(1,:), ...
    velocity_points(2,:),50,'filled');
axis(ax2,'equal'); grid(ax2,'on'); box(ax2,'on');
xlim(ax2,domain_scale*[-1 1]); ylim(ax2,domain_scale*[-1 1]);
xlabel(ax2,'dq_1'); ylabel(ax2,'dq_2');
title(ax2,sprintf('Representative velocity samples (%d)', ...
    velocity_count), ...
    'FontWeight','normal');
title(layout,sprintf('Structured inducing set: %d positions x %d velocities = %d states', ...
    position_count,velocity_count,size(InducingPoints_Coordinates,2)), ...
    'FontWeight','bold');
exportgraphics(fig,strrep(output_file,'.mat','.png'),'Resolution',220);
close(fig);

fprintf('Saved structured inducing points: %s\n',output_file);
end
