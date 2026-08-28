function [TrackingError_vector, t_set] = run_simulation_inducing_point( ...
    CurrentMode, SaveFolderName, SaveFileName, use_formation, ...
    UnknownScaleOverride, DisturbanceScaleOverride, OfflineDataQuantityOverride, ...
    InitialSeedOverride)
% Pure inducing-point online-learning simulation.
%
% Purpose of this version:
%   Keep only the core algorithmic logic for advisor discussion.
%
% Core logic:
%   1) Offline training data = 0.
%   2) Every 0.1 s, every agent adds its own current state/dynamics pair to
%      its own LocalGP.  Thus a 30 s run gives 300 online samples per agent,
%      matching the shared code's per-agent dataset ownership/counting.
%      The 0.1 s acquisition itself is time-triggered in this experiment.
%   3) After each online update, all LocalGPs are projected onto a fixed
%      inducing-point set P.
%   4) IP-DAC remains active for the full simulation. At every physical
%      step, equation (17) decides broadcasts and a fixed communication
%      budget is used. Zeta and the last broadcasts persist across steps.
%      There is no tracking/consensus terminal condition.
%   5) IP-AC restarts when the inducing-point information P is refreshed.
%      Xi^(0)=P_current, while the last values actually broadcast to
%      neighbors persist across windows. The fixed R periodic ET-AC rounds
%      are spread over the following real 0.1 s physical-time window; the
%      controller therefore uses the current intermediate Xi at each step.
%   5) A MaskedGP is rebuilt from the current (not terminal) DAC state.
%   6) For inducing-point modes, the controller always uses this shared
%      inducing-point MaskedGP prediction. There is no shadow/local switch.
%
% Supported modes:
%   poe, gpoe, moe, bcm, rbcm                 : IP-DAC
%   poe_offline, gpoe_offline, moe_offline,
%   bcm_offline, rbcm_offline                 : offline IP-DAC
%   poe_ac, gpoe_ac, moe_ac, bcm_ac, rbcm_ac : IP-AC
%   poe_ac_offline, ..., rbcm_ac_offline      : offline IP-AC
%   local                                    : LocalGP only baseline
%   without_gp                               : zero compensation baseline
%   without_gp_unknown0                      : nominal baseline with no true unknown term
%   exact                                    : exact unknown dynamics baseline


if nargin < 4
    use_formation = true;
end
if nargin < 5 || isempty(UnknownScaleOverride)
    UnknownScaleOverride = [];
end
if nargin < 6 || isempty(DisturbanceScaleOverride)
    DisturbanceScaleOverride = [];
end
if nargin < 7 || isempty(OfflineDataQuantityOverride)
    OfflineDataQuantityOverride = [];
end
if nargin < 8 || isempty(InitialSeedOverride)
    InitialSeedOverride = [];
end
if nargin < 2 || isempty(SaveFolderName)
    SaveFolderName = '';
end
if nargin < 3 || isempty(SaveFileName)
    SaveFileName = '';
end
if ~isempty(SaveFolderName) && ~isfolder(SaveFolderName)
    mkdir(SaveFolderName);
end

rng(0);

%% Fixed algorithm parameters
ACFixedIterations = get_env_positive_integer( ...
    'CONTROL_AC_FIXED_ITERATIONS',10);
% Kept in saved result files for backward-compatible diagnostics.  Online
% learning is deliberately fixed to the per-agent dataset convention.
OnlineAgentPolicy = 'all_agents';

UnknownScale = 0.2;
% DisturbanceScale scales only the GP observation noise.
DisturbanceScale = 0.1;
if ~isempty(UnknownScaleOverride)
    UnknownScale = UnknownScaleOverride;
end
if ~isempty(DisturbanceScaleOverride)
    DisturbanceScale = DisturbanceScaleOverride;
end

%% 1. System parameters
SystemOrder = 2;
q_dim = 2;
x_dim = q_dim * SystemOrder;
y_dim = q_dim;

m1 = 1;
m2 = 1;
L1 = 1;
L2 = 1;
AgentQuantity = 6;
LeaderQuantity = 1;

%% 2. Topology
MultiAgentSystem = Manipulator_2D_2DoF_SetMASTopology(AgentQuantity, LeaderQuantity);
L = MultiAgentSystem.Agent_Topology.LaplacianMatrix;
B = MultiAgentSystem.Agent_Leader_Topology.ConnectionMatrix(:,1);
L_lap = L;
ACPeriodicSigma = get_env_positive_scalar( ...
    'CONTROL_AC_PERIODIC_SIGMA',0.3);
ACForceFullBroadcast = strcmp(strtrim(getenv( ...
    'CONTROL_AC_FORCE_FULL_BROADCAST')),'1');
% Diagnostics store the two scalar sides of the paper's squared,
% agent-level detector directly:
%   ||Xi_i-XiHat_i||_F^2/(pM)
%   versus sigma_i^2*sum_j a_ij||XiHat_i-XiHat_j||_F^2/(pM).
ACDetectorStoresSquared = true;
DACTriggerEpsilon = get_env_positive_scalar( ...
    'CONTROL_DAC_TRIGGER_EPSILON', 5e-3);

%% 3. Controller parameters
c = 10;
lambda_set = [1; 1];
lambda_n = lambda_set(end);
lambda_vector = lambda_set(1:SystemOrder-1);

Lambda = [zeros(SystemOrder-2,1), eye(SystemOrder-2); ...
          -lambda_set(1)/lambda_n, -lambda_set(2:end-1)/lambda_n];

Qes = eye(SystemOrder-1);
Pes = care(Lambda, [], Qes);
Pe  = kron(Pes, eye(AgentQuantity));
Qe  = kron(Qes, eye(AgentQuantity)); %#ok<NASGU>

q  = (L + diag(B)) \ ones(AgentQuantity,1);
Pr = diag(1./q);
Qr = Pr*(L+diag(B)) + (L+diag(B))'*Pr;

t_vec = [zeros(SystemOrder-2,1); 1/lambda_n];
Phi = Pr * kron(lambda_vector' * t_vec, eye(AgentQuantity));
Psi = Pr * kron(lambda_vector' * Lambda, eye(AgentQuantity)) + ...
      kron(t_vec' * Pes, eye(AgentQuantity));
Qz = [c*lambda_n*Qr - 2*Phi, -Psi; -Psi', Pe];

if all(real(eig(Qz)) > 0) && all(real(eig(Lambda)) < 0)
    fprintf('The controller is stable!\n');
else
    error('Controller is not stable.');
end

%% 4. Time

t_start = 0;
t_end = 10;
t_end = get_env_positive_scalar('CONTROL_SIM_END_TIME', t_end);
t_step = 0.01;
t_set = t_start:t_step:t_end;
T = numel(t_set);
IPOnlineUpdateInterval = 0.1;
IPOnlineUpdateStep = max(1, round(IPOnlineUpdateInterval / t_step));
ACOnlineUpdateInterval = IPOnlineUpdateInterval;
ACOnlineUpdateStep = IPOnlineUpdateStep;
EnableProjectionDiagnostics = strcmp(getenv('CONTROL_IP_PROJECTION_DIAGNOSTICS'), '1');
EnableACTriggerDiagnostics = strcmp( ...
    getenv('CONTROL_AC_TRIGGER_DIAGNOSTICS'),'1');
AggregationParameters = control_aggregation_parameters();
RBCMBetaMax = AggregationParameters.rbcm_beta_max;

%% 5. Reference trajectory and formation offsets
[xl_set, xlr_set, ~] = Manipulator_2D_2DoF_LeaderDynamics(t_set, L1);

s_all_set  = nan(x_dim*AgentQuantity, T);
sr_all_set = nan(q_dim*AgentQuantity, T);
for AgentNr = 1:AgentQuantity
    [s_all_set((AgentNr-1)*x_dim+(1:x_dim),:), ...
     sr_all_set((AgentNr-1)*q_dim+(1:q_dim),:)] = ...
        Manipulator_2D_2DoF_RelativePositionDynamics(t_set, AgentNr, AgentQuantity);
end

if ~use_formation
    s_all_set  = zeros(size(s_all_set));
    sr_all_set = zeros(size(sr_all_set));
end

%% 6. LocalGPs: zero offline data by default
SigmaF = 1;
% Separate spatial and velocity length scales.  Defaults reproduce the
% previous isotropic kernel; environment overrides are useful when testing
% a position-dense / velocity-sparse inducing-point design.
PositionLengthScale = get_env_positive_scalar( ...
    'CONTROL_IP_POSITION_LENGTH_SCALE',0.5);
VelocityLengthScale = get_env_positive_scalar( ...
    'CONTROL_IP_VELOCITY_LENGTH_SCALE',0.5);
SigmaL = [PositionLengthScale*ones(2,1); ...
          VelocityLengthScale*ones(2,1)];
SigmaN = 0.01;
% Numerical jitter used only when reconstructing a MaskedGP from the
% inducing-point targets. It is separate from the local-GP noise above.
IPReconstructionSigmaN = get_env_positive_scalar( ...
    'CONTROL_IP_RECONSTRUCTION_SIGMA_N',1e-6);
GP_tau = 1e-8;
GP_delta = 0.1;
DomainScale = 1.5;
X_min = DomainScale * [-1,-1,-1,-1];
X_max = DomainScale * [ 1, 1, 1, 1];

MaxDataQuantity = 600;
DefaultOfflineDataQuantity = 350;

mode_lower = lower(CurrentMode);
is_without_gp_unknown0_mode = strcmpi(mode_lower, 'without_gp_unknown0');
if is_without_gp_unknown0_mode
    UnknownScale = 0;
end
dac_methods = {'poe','gpoe','moe','bcm','rbcm'};
ac_methods  = {'poe_ac','gpoe_ac','moe_ac','bcm_ac','rbcm_ac'};
offline_methods = strcat(dac_methods, '_offline');
offline_ac_methods = strcat(ac_methods, '_offline');
is_offline_ip_mode = ismember(mode_lower, [offline_methods, offline_ac_methods]);
is_dac_mode = ismember(mode_lower, [dac_methods, offline_methods]);
is_ac_mode  = ismember(mode_lower, [ac_methods, offline_ac_methods]);
if is_offline_ip_mode
    if isempty(OfflineDataQuantityOverride)
        OfflineDataQuantity = DefaultOfflineDataQuantity;
    else
        OfflineDataQuantity = OfflineDataQuantityOverride;
    end
else
    OfflineDataQuantity = 0;
end

LocalGP_set = cell(AgentQuantity, 1);
for AgentNr = 1:AgentQuantity
    LocalGP_set{AgentNr} = LocalGP_MultiOutput(x_dim, y_dim, MaxDataQuantity, ...
        SigmaN, SigmaF, SigmaL);
    LocalGP_set{AgentNr}.tau   = GP_tau;
    LocalGP_set{AgentNr}.delta = GP_delta;
    LocalGP_set{AgentNr}.xMax  = X_max;
    LocalGP_set{AgentNr}.xMin  = X_min;
end

if is_offline_ip_mode
    for AgentNr = 1:AgentQuantity
        X_in = 2 * (rand(x_dim, OfflineDataQuantity) - 0.5) * DomainScale;
        Y_in = UnknownScale * Manipulator_2D_2DoF_UnknownDynamics(X_in) + ...
            DisturbanceScale * LocalGP_set{AgentNr}.SigmaN * ...
            randn(y_dim, OfflineDataQuantity);
        LocalGP_set{AgentNr}.add_Alldata(X_in, Y_in);
    end
end

if is_offline_ip_mode
    fprintf('Offline data: %d per agent. Online update: disabled.\n', OfflineDataQuantity);
elseif is_ac_mode
    fprintf('Offline data: %d per agent. IP online update interval: %.3f s.\n', ...
        OfflineDataQuantity, IPOnlineUpdateInterval);
else
    fprintf('Offline data: %d per agent. IP online update interval: %.3f s.\n', ...
        OfflineDataQuantity, IPOnlineUpdateInterval);
end
fprintf('UnknownScale = %.3f, DisturbanceScale = %.3f.\n', UnknownScale, DisturbanceScale);

%% 7. Mode and inducing-point setup
is_ip_mode  = is_dac_mode || is_ac_mode;
is_local_mode = strcmpi(mode_lower, 'local');
is_exact_mode = strcmpi(mode_lower, 'exact');
is_without_gp_mode = strcmpi(mode_lower, 'without_gp') || is_without_gp_unknown0_mode;

if is_dac_mode
    base_method = strrep(mode_lower, '_offline', '');
    ConsensusInputTriggerTol = NaN;
elseif is_ac_mode
    base_method = strrep(strrep(mode_lower, '_offline', ''), '_ac', '');
    ConsensusInputTriggerTol = NaN;
elseif is_local_mode || is_exact_mode || is_without_gp_mode
    ConsensusInputTriggerTol = NaN;
    if is_without_gp_mode
        base_method = 'without_gp';
    else
        base_method = mode_lower;
    end
else
    error('Unknown mode: %s', CurrentMode);
end

% Keep DAC unchanged. IP-AC performs a fixed number of communication
% checks over each real 0.1 s P-reference window. The checks are no longer
% executed instantaneously when P changes: they are distributed over the
% physical integration steps in that window.
Kappa_P_DAC = 1;
LegacyACConsensusStep = get_env_positive_scalar( ...
    'CONTROL_AC_CONSENSUS_STEP',0.2);
ACChecksPerPhysicalStepExact = ACFixedIterations/ACOnlineUpdateStep;
ACChecksPerPhysicalStep = round(ACChecksPerPhysicalStepExact);
if ACChecksPerPhysicalStep < 1 || ...
        abs(ACChecksPerPhysicalStep-ACChecksPerPhysicalStepExact) > 1e-12
    error(['ACFixedIterations must be an integer multiple of the ten ' ...
        'physical steps in one 0.1 s P-update window.']);
end
ACCommunicationInterval = t_step/ACChecksPerPhysicalStep;
Kappa_P_AC = get_env_positive_scalar( ...
    'CONTROL_AC_KAPPA_P',LegacyACConsensusStep/t_step);
ACConsensusStep = ACCommunicationInterval*Kappa_P_AC;
ACTimingMode = 'physical-time-fixed-rounds';
NumInducingPoints = 400;
num_inducing_override = str2double(getenv('CONTROL_IP_NUM_INDUCING_POINTS'));
if isfinite(num_inducing_override) && num_inducing_override > 0
    NumInducingPoints = round(num_inducing_override);
end
inducing_point_file = strtrim(getenv('CONTROL_IP_INDUCING_POINT_FILE'));
if isempty(inducing_point_file)
    % Generate Z from its own deterministic stream.  This keeps the same
    % inducing-point locations when only AC/DAC parameters are swept, and
    % does not alter the random stream used later for the initial state.
    InducingPointSeed = get_env_positive_integer( ...
        'CONTROL_IP_INDUCING_POINT_SEED',42);
    rng_state_before_inducing_points = rng;
    rng(InducingPointSeed,'twister');
    InducingPoints_Coordinates = ...
        2*DomainScale*rand(x_dim, NumInducingPoints) - DomainScale;
    rng(rng_state_before_inducing_points);
    InducingPointSource = sprintf( ...
        'uniform-domain-seed-%d',InducingPointSeed);
else
    InducingPointSeed = nan;
    inducing_point_data = load(inducing_point_file,'InducingPoints_Coordinates');
    if ~isfield(inducing_point_data,'InducingPoints_Coordinates')
        error('Inducing-point file must contain InducingPoints_Coordinates.');
    end
    InducingPoints_Coordinates = inducing_point_data.InducingPoints_Coordinates;
    if size(InducingPoints_Coordinates,1) ~= x_dim
        error('Loaded inducing points must have %d rows.',x_dim);
    end
    NumInducingPoints = size(InducingPoints_Coordinates,2);
    if NumInducingPoints < 1 || any(~isfinite(InducingPoints_Coordinates(:)))
        error('Loaded inducing-point set is empty or contains nonfinite values.');
    end
    InducingPointSource = inducing_point_file;
end
fprintf('Inducing points: M=%d, source=%s.\n', ...
    NumInducingPoints,InducingPointSource);

DACStepsPerPhysicalStep = get_env_positive_integer( ...
    'CONTROL_DAC_STEPS_PER_PHYSICAL_STEP',1);

%% 8. Initial state
if ~isempty(InitialSeedOverride)
    rng(InitialSeedOverride);
elseif ~is_offline_ip_mode
    rng(42);
end
x_all = rand(x_dim*AgentQuantity, 1);
x_all_set = nan(x_dim*AgentQuantity, T);
x_all_set(:,1) = x_all;

vartheta_all_set = nan(x_dim*AgentQuantity, T);
vartheta_all_set(:,1) = x_all - s_all_set(:,1) - kron(ones(AgentQuantity,1), xl_set(:,1));

TrackingError_vector = zeros(1,T);
f_hat_matrix = zeros(y_dim, AgentQuantity);
f_true_matrix = zeros(y_dim, AgentQuantity);
f_direct_matrix = nan(y_dim, AgentQuantity);
f_ip_cen_matrix = nan(y_dim, AgentQuantity);

prediction_error_norm_vector = nan(1,T);
f_hat_norm_vector = nan(1,T);
direct_prediction_error_norm_vector = nan(1,T);
ip_projection_gap_norm_vector = nan(1,T);
f_hat_all_set = nan(y_dim, AgentQuantity, T);
f_true_all_set = nan(y_dim, AgentQuantity, T);
f_direct_all_set = nan(y_dim, AgentQuantity, T);
f_ip_cen_all_set = nan(y_dim, AgentQuantity, T);
ip_cen_prediction_error_norm_vector = nan(1,T);
aggregation_min_abs_den_vector = nan(1,T);
aggregation_min_signed_den_vector = nan(1,T);
aggregation_negative_den_count_vector = nan(1,T);
aggregation_small_den_count_vector = nan(1,T);
aggregation_max_abs_phi_vector = nan(1,T);
online_update_count = zeros(AgentQuantity,1);
online_learning_trigger_count_set = zeros(AgentQuantity,T-1);
online_learning_error_set = nan(AgentQuantity,T-1);
projection_update_set = zeros(1,T);
consensus_communication_update_set = zeros(1,T);
consensus_input_change_vector = nan(1,T);
consensus_input_change_per_agent_set = nan(AgentQuantity,T);
consensus_input_numerator_change_per_agent_set = nan(AgentQuantity,T);
consensus_input_denominator_change_per_agent_set = nan(AgentQuantity,T);
projection_step_change_per_agent_set = nan(AgentQuantity,T);
projection_absolute_change_rms_set = nan(1,T);
projection_relative_change_per_agent_set = nan(AgentQuantity,T);
projection_step_numerator_change_per_agent_set = nan(AgentQuantity,T);
projection_step_denominator_change_per_agent_set = nan(AgentQuantity,T);
consensus_input_trigger_count_set = zeros(AgentQuantity,T-1);

MaskedGP = [];
CenMaskedGP = [];
Xi_final = [];
P_inducing = [];
P_consensus_reference = [];
P_previous_projection = [];
p_dim = [];
Zeta_dac = [];
Xi_last_trigger = [];
Xi_hat_ac = [];
Xi_ac = [];
P_average_ac = [];
ac_round_in_window = 0;
ac_solve_count = 0;
dac_event_count_per_agent = zeros(AgentQuantity,1);
dac_broadcast_event_count_per_agent = zeros(AgentQuantity,1);
ac_event_count_per_agent = zeros(AgentQuantity,1);
dac_broadcast_count_set = zeros(AgentQuantity,T-1);
dac_trigger_count_per_agent_point_set = zeros( ...
    AgentQuantity,NumInducingPoints,T-1,'uint16');
dac_inner_broadcast_count_set = zeros(AgentQuantity,T-1);
dac_iteration_count_set = zeros(1,T-1);
dac_tracking_error_set = nan(1,T-1);
dac_tracking_error_per_agent_set = nan(AgentQuantity,T-1);
dac_consensus_disagreement_set = nan(1,T-1);
dac_inner_trigger_instance_set = false( ...
    AgentQuantity,DACStepsPerPhysicalStep,T-1);
ac_broadcast_count_set = zeros(AgentQuantity,T-1);
ac_inner_broadcast_count_set = zeros(AgentQuantity,T-1);
ac_physical_trigger_count_set = zeros(AgentQuantity,T-1);
ac_iteration_count_set = zeros(1,T-1);
ac_consensus_disagreement_set = nan(1,T-1);
ac_inner_trigger_instance_set = false( ...
    AgentQuantity,ACChecksPerPhysicalStep,T-1);
ac_tracking_error_history_set = cell(1,T-1);
if EnableACTriggerDiagnostics
    ac_trigger_measure_set = nan( ...
        AgentQuantity,ACChecksPerPhysicalStep,T-1,'single');
    ac_trigger_threshold_set = nan( ...
        AgentQuantity,ACChecksPerPhysicalStep,T-1,'single');
else
    ac_trigger_measure_set = single([]);
    ac_trigger_threshold_set = single([]);
end

if is_ip_mode && ~is_offline_ip_mode
    fprintf('Mode: %s. Projection refresh: every %.3f s.\n', ...
        CurrentMode, IPOnlineUpdateInterval);
elseif is_ip_mode
    fprintf('Mode: %s. Projection refresh: offline initialization only.\n', ...
        CurrentMode);
else
    fprintf('Mode: %s. Projection refresh: every simulation step.\n', CurrentMode);
end
if is_ip_mode
    if is_dac_mode
        if is_offline_ip_mode
            fprintf('Control source: offline LocalGPs aggregated by IP-DAC.\n');
        else
            fprintf('Control source: per-agent MaskedGP driven by IP-DAC.\n');
        end
    else
        fprintf('Control source: shared inducing-point MaskedGP.\n');
    end
elseif is_local_mode
    fprintf('Control source: LocalGP baseline.\n');
elseif is_without_gp_mode
    if is_without_gp_unknown0_mode
        fprintf('Control source: zero compensation nominal baseline; true unknown dynamics are disabled.\n');
    else
        fprintf('Control source: zero compensation baseline; true unknown dynamics remain in the plant.\n');
    end
elseif is_exact_mode
    fprintf('Control source: exact unknown dynamics.\n');
end

%% 9. Main loop
opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-3, ...
    'InitialStep', t_step, 'MaxStep', t_step, 'Refine', 1);
elapsed_tic = tic;

for t_Nr = 1:T-1
    t = t_set(t_Nr);

    x_l_r = xlr_set(:, t_Nr);
    x_all = x_all_set(:, t_Nr);
    x_all_matrix = reshape(x_all, x_dim, AgentQuantity);
    x_all_cell = ET_MAS_GP_Leader_vector2cell(x_all, AgentQuantity, 1);

    s_all = s_all_set(:, t_Nr);
    s_r_all = sr_all_set(:, t_Nr);
    s_r_cell = ET_MAS_GP_Leader_vector2cell(s_r_all, AgentQuantity, 1);

    x_tilde_all = x_all - s_all;
    x_tilde_cell = ET_MAS_GP_Leader_vector2cell(x_tilde_all, AgentQuantity, SystemOrder);

    vartheta_all = vartheta_all_set(:, t_Nr);
    vartheta_cell = ET_MAS_GP_Leader_vector2cell(vartheta_all, AgentQuantity, SystemOrder);
    TrackingError_vector(t_Nr) = norm(vartheta_all);
    is_ip_online_step = is_ip_mode && (mod(t_Nr-1, IPOnlineUpdateStep) == 0);
    is_learning_update_step = (is_ip_mode && is_ip_online_step) || is_local_mode;
    step_online_trigger_count = zeros(AgentQuantity,1);
    has_online_learning_update = false;

    %% 9.1 Per-agent time-triggered LocalGP online update
    if is_ip_online_step && ~is_offline_ip_mode
        % Each agent owns a separate GP dataset.  At every 0.1 s learning
        % instant, all agents add one local data pair; there is no global
        % round-robin scheduler and no competition for a network-wide slot.
        for AgentNr = 1:AgentQuantity
            agent_state_idx = (AgentNr-1)*x_dim + (1:x_dim);
            agent_tracking_error = norm(vartheta_all(agent_state_idx));
            online_learning_error_set(AgentNr,t_Nr) = agent_tracking_error;

            x_i = x_all_matrix(:,AgentNr);
            y_i = UnknownScale * Manipulator_2D_2DoF_UnknownDynamics(x_i) + ...
                  DisturbanceScale * LocalGP_set{AgentNr}.SigmaN * randn(y_dim,1);

            if LocalGP_set{AgentNr}.DataQuantity >= LocalGP_set{AgentNr}.MaxDataQuantity
                LocalGP_set{AgentNr}.downdateParam(1);
            end
            LocalGP_set{AgentNr}.addPoint(x_i, y_i);
            online_update_count(AgentNr) = online_update_count(AgentNr) + 1;
            step_online_trigger_count(AgentNr) = 1;
            has_online_learning_update = true;
        end
        online_learning_trigger_count_set(:,t_Nr) = step_online_trigger_count;
    elseif is_learning_update_step && ~is_offline_ip_mode
        for AgentNr = 1:AgentQuantity
            x_i = x_all_matrix(:,AgentNr);
            y_i = UnknownScale * Manipulator_2D_2DoF_UnknownDynamics(x_i) + ...
                  DisturbanceScale * LocalGP_set{AgentNr}.SigmaN * randn(y_dim,1);

            if LocalGP_set{AgentNr}.DataQuantity >= LocalGP_set{AgentNr}.MaxDataQuantity
                LocalGP_set{AgentNr}.downdateParam(1);
            end
            LocalGP_set{AgentNr}.addPoint(x_i, y_i);
            online_update_count(AgentNr) = online_update_count(AgentNr) + 1;
        end
    end

    %% 9.2 Inducing-point projection and ET consensus
    if is_ip_mode
        should_refresh_projection = (~is_offline_ip_mode && is_ip_online_step && ...
            has_online_learning_update) || ...
            isempty(P_inducing);
        if should_refresh_projection
            projection_update_set(t_Nr) = 1;

            [P_inducing, p_dim, var_diag] = gp_masked_aggregation_init( ...
                LocalGP_set, AgentQuantity, NumInducingPoints, ...
                InducingPoints_Coordinates, base_method);

            % Diagnostic: track how close each agent's local GP posterior
            % variance (at the inducing points) is getting to the
            % posterior_var_floor over time. If min_var keeps dropping
            % toward the floor, the precision channel (~1/var) it feeds
            % into P can jump sharply -- this is the suspected driver of
            % the large P[k] jumps seen around specific times.
            if ~exist('var_diag_min_set','var')
                var_diag_min_set = nan(AgentQuantity, T);
                var_diag_near_floor_set = nan(AgentQuantity, T);
            end
            var_diag_min_set(:,t_Nr) = var_diag.min_var;
            var_diag_near_floor_set(:,t_Nr) = var_diag.near_floor_count;

            % Shadow IP-CEN reference: use the exact average of the current
            % information vectors, without DAC, AC, or event triggering.
            % It is evaluated on the same closed-loop states as IP-DAC and
            % never enters the controller.
            if EnableProjectionDiagnostics
                Xi_ip_cen = mean(P_inducing,2);
                CenMaskedGP = build_or_update_shared_maskedgp_from_xi( ...
                    Xi_ip_cen, base_method, AgentQuantity, ...
                    NumInducingPoints, p_dim, InducingPoints_Coordinates, ...
                    SigmaF, SigmaL, x_dim, y_dim, CenMaskedGP);
            end

            if ~isempty(P_previous_projection)
                P_step_delta = P_inducing - P_previous_projection;
                projection_absolute_change_rms_set(t_Nr) = ...
                    sqrt(mean(P_step_delta(:).^2));
                projection_step_change_per_agent_set(:,t_Nr) = ...
                    rms_change_per_agent(P_step_delta, 1:p_dim);
                previous_P_rms = ...
                    rms_change_per_agent(P_previous_projection,1:p_dim);
                projection_relative_change_per_agent_set(:,t_Nr) = ...
                    projection_step_change_per_agent_set(:,t_Nr) ./ ...
                    max(previous_P_rms,eps);
                projection_step_numerator_change_per_agent_set(:,t_Nr) = ...
                    rms_change_per_agent(P_step_delta, 1:2:p_dim);
                projection_step_denominator_change_per_agent_set(:,t_Nr) = ...
                    rms_change_per_agent(P_step_delta, 2:2:p_dim);
            end
            P_previous_projection = P_inducing;
        end

        if is_dac_mode
            % Zeta and Xi_last_trigger are persistent DAC states.  The DAC
            % never declares convergence and never resets when P changes.
            % Instead, one fixed DAC step is taken at every physical step;
            % equation (17) alone decides whether communication occurs.
            is_first_dac_step = isempty(Zeta_dac);
            if is_first_dac_step
                Zeta_dac = zeros(size(P_inducing));
                Xi_last_trigger = P_inducing;
            end

            step_trigger_count_per_agent_point = zeros( ...
                AgentQuantity,NumInducingPoints,'uint16');

            % DACStepsPerPhysicalStep is a numerical subdivision of one
            % physical interval.  Keep the total DAC evolution over that
            % interval equal to t_step when comparing different round
            % counts; otherwise R rounds would accelerate DAC by R times.
            dac_step_size = t_step/DACStepsPerPhysicalStep;
            for dac_step = 1:DACStepsPerPhysicalStep
                L_Xi_hat = laplacian_multiply_agent_dim_local( ...
                    Xi_last_trigger,L_lap);
                Zeta_dac = Zeta_dac+dac_step_size*Kappa_P_DAC*L_Xi_hat;
                Xi_now = P_inducing-Zeta_dac;

                force_initial_broadcast = ...
                    is_first_dac_step && dac_step == 1;
                [Xi_last_trigger, triggered_points] = ...
                    dac_equation17_broadcast_update( ...
                    Xi_now, Xi_last_trigger, L_lap, ...
                    DACTriggerEpsilon, force_initial_broadcast);

                step_trigger_count_per_agent_point = ...
                    step_trigger_count_per_agent_point ...
                    + uint16(triggered_points);
                dac_inner_trigger_instance_set(:,dac_step,t_Nr) = ...
                    any(triggered_points,2);
            end

            % These errors are diagnostics only.  Neither can stop DAC.
            Xi_final = Xi_now;
            P_average = mean(P_inducing,2);
            tracking_difference = Xi_now-P_average;
            tracking_error_norm = sqrt(sum(tracking_difference.^2,1));
            dac_tracking_error_set(t_Nr) = max(tracking_error_norm(:));
            dac_tracking_error_per_agent_set(:,t_Nr) = reshape( ...
                max(tracking_error_norm,[],3),AgentQuantity,1);

            consensus_disagreement = 0;
            for agent_i = 1:AgentQuantity
                neighbor_list = find(L_lap(agent_i,:) < 0);
                for agent_j = neighbor_list
                    Xi_difference = Xi_now(:,agent_i,:) ...
                        - Xi_now(:,agent_j,:);
                    Xi_difference_norm = sqrt(sum(Xi_difference.^2,1));
                    consensus_disagreement = max(consensus_disagreement, ...
                        max(Xi_difference_norm(:)));
                end
            end
            dac_consensus_disagreement_set(t_Nr) = consensus_disagreement;
            dac_iteration_count_set(t_Nr) = DACStepsPerPhysicalStep;

            step_trigger_count = mean( ...
                double(step_trigger_count_per_agent_point),2);
            step_broadcast_count = sum( ...
                double(dac_inner_trigger_instance_set(:,:,t_Nr)),2);
            dac_trigger_count_per_agent_point_set(:,:,t_Nr) = ...
                step_trigger_count_per_agent_point;
            dac_event_count_per_agent = dac_event_count_per_agent ...
                + sum(double(step_trigger_count_per_agent_point),2);
            dac_broadcast_event_count_per_agent = ...
                dac_broadcast_event_count_per_agent + step_broadcast_count;
            dac_broadcast_count_set(:,t_Nr) = step_broadcast_count;
            dac_inner_broadcast_count_set(:,t_Nr) = step_trigger_count;
            consensus_communication_update_set(t_Nr) = ...
                any(step_trigger_count_per_agent_point(:) > 0);

            MaskedGP = build_agent_maskedgp_from_xi(Xi_final, base_method, ...
                AgentQuantity, NumInducingPoints, p_dim, InducingPoints_Coordinates, ...
                SigmaF, SigmaL, x_dim, y_dim, MaskedGP);

        else
            %% IP-AC: fixed R checks spread over the real 0.1 s window
            % A new P reference arrives every 0.1 s. Xi is reset to that P,
            % then exactly ACFixedIterations ET/AC rounds are executed over
            % the following physical steps. There is no tolerance loop and
            % no extra communication after the fixed round budget.
            step_input_trigger_count = zeros(AgentQuantity,1);
            if should_refresh_projection && ~isempty(P_consensus_reference)
                P_delta = P_inducing-P_consensus_reference;
                consensus_input_change_vector(t_Nr) = ...
                    sqrt(mean(P_delta(:).^2,'omitnan'));
                for AgentNr = 1:AgentQuantity
                    P_delta_i = P_delta(:,AgentNr,:);
                    consensus_input_change_per_agent_set(AgentNr,t_Nr) = ...
                        sqrt(mean(P_delta_i(:).^2,'omitnan'));
                end
                consensus_input_numerator_change_per_agent_set(:,t_Nr) = ...
                    rms_change_per_agent(P_delta,1:2:p_dim);
                consensus_input_denominator_change_per_agent_set(:,t_Nr) = ...
                    rms_change_per_agent(P_delta,2:2:p_dim);
                step_input_trigger_count(:) = 1;
            elseif should_refresh_projection
                step_input_trigger_count(:) = 1;
            end

            if should_refresh_projection
                consensus_input_trigger_count_set(:,t_Nr) = ...
                    step_input_trigger_count;
                P_consensus_reference = P_inducing;
                Xi_ac = P_consensus_reference;
                P_average_ac = mean(P_consensus_reference,2);
                ac_round_in_window = 0;
                ac_solve_count = ac_solve_count+1;
                if isempty(Xi_hat_ac)
                    Xi_hat_ac = Xi_ac;
                end
            end

            if ~isempty(P_consensus_reference) && ...
                    ac_round_in_window < ACFixedIterations
                checks_this_step = min(ACChecksPerPhysicalStep, ...
                    ACFixedIterations-ac_round_in_window);
                step_trigger_count = zeros(AgentQuantity,1);
                step_trigger_measure_history = nan( ...
                    AgentQuantity,ACChecksPerPhysicalStep,'single');
                step_trigger_threshold_history = nan( ...
                    AgentQuantity,ACChecksPerPhysicalStep,'single');
                tracking_error_history = nan(1,checks_this_step+1);
                initial_average_error = Xi_ac-P_average_ac;
                tracking_error_history(1) = ...
                    sqrt(mean(initial_average_error(:).^2));

                for check_i = 1:checks_this_step
                    ac_round_in_window = ac_round_in_window+1;

                    % Agent-level PETC detector, evaluated at the real
                    % communication time separated by ACCommunicationInterval.
                    ac_lhs = zeros(AgentQuantity,1);
                    ac_rhs = zeros(AgentQuantity,1);
                    ac_trigger_mask = false(AgentQuantity,1);
                    for agent_i = 1:AgentQuantity
                        neighbor_list = find(L_lap(agent_i,:) < 0);
                        local_error_i = Xi_ac(:,agent_i,:) ...
                            - Xi_hat_ac(:,agent_i,:);
                        packet_entry_count = numel(local_error_i);
                        ac_lhs_i = sum(local_error_i.^2,'all') ...
                            /packet_entry_count;

                        neighbor_term_i = 0;
                        for agent_j = neighbor_list
                            a_ij = -L_lap(agent_i,agent_j);
                            cached_gap_ij = Xi_hat_ac(:,agent_i,:) ...
                                - Xi_hat_ac(:,agent_j,:);
                            neighbor_term_i = neighbor_term_i ...
                                + a_ij*sum(cached_gap_ij.^2,'all') ...
                                /packet_entry_count;
                        end

                        ac_rhs_i = ACPeriodicSigma^2*neighbor_term_i;
                        ac_lhs(agent_i) = ac_lhs_i;
                        ac_rhs(agent_i) = ac_rhs_i;
                        ac_trigger_mask(agent_i) = ...
                            ACForceFullBroadcast ...
                            || (ac_solve_count == 1 && ...
                            ac_round_in_window == 1) ...
                            || ac_lhs_i >= ac_rhs_i;
                    end

                    % All broadcasts at this check are simultaneous.
                    for agent_i = find(ac_trigger_mask).'
                        Xi_hat_ac(:,agent_i,:) = Xi_ac(:,agent_i,:);
                    end
                    step_trigger_count = step_trigger_count ...
                        + double(ac_trigger_mask);
                    ac_inner_trigger_instance_set(:,check_i,t_Nr) = ...
                        ac_trigger_mask;
                    step_trigger_measure_history(:,check_i) = single(ac_lhs);
                    step_trigger_threshold_history(:,check_i) = single(ac_rhs);

                    % One Euler AC update over the true communication interval.
                    L_Xi_hat = laplacian_multiply_agent_dim_local( ...
                        Xi_hat_ac,L_lap);
                    Xi_ac = Xi_ac-ACCommunicationInterval*Kappa_P_AC*L_Xi_hat;
                    average_error = Xi_ac-P_average_ac;
                    tracking_error_history(check_i+1) = ...
                        sqrt(mean(average_error(:).^2));
                end

                ac_iteration_count_set(t_Nr) = checks_this_step;
                ac_disagreement = abs(Xi_ac-P_average_ac);
                ac_consensus_disagreement_set(t_Nr) = ...
                    mean(max(ac_disagreement,[],2),'all');
                ac_tracking_error_history_set{t_Nr} = tracking_error_history;
                if EnableACTriggerDiagnostics
                    ac_trigger_measure_set(:,:,t_Nr) = ...
                        step_trigger_measure_history;
                    ac_trigger_threshold_set(:,:,t_Nr) = ...
                        step_trigger_threshold_history;
                end

                Xi_final = Xi_ac;
                step_physical_trigger_count = double(step_trigger_count > 0);
                ac_event_count_per_agent = ac_event_count_per_agent ...
                    + step_trigger_count;
                ac_broadcast_count_set(:,t_Nr) = step_trigger_count;
                ac_inner_broadcast_count_set(:,t_Nr) = step_trigger_count;
                ac_physical_trigger_count_set(:,t_Nr) = ...
                    step_physical_trigger_count;
                consensus_communication_update_set(t_Nr) = ...
                    any(step_trigger_count > 0);
            end

            % The controller sees the current, generally not-yet-converged
            % AC state at every physical step. Reuse the fixed kernel
            % factorization and update only the inducing targets.
            if ~isempty(Xi_ac)
                Xi_final = Xi_ac;
                MaskedGP = build_or_update_shared_maskedgp_from_xi( ...
                    Xi_final,base_method,AgentQuantity,NumInducingPoints, ...
                    p_dim,InducingPoints_Coordinates,SigmaF,SigmaL, ...
                    x_dim,y_dim,MaskedGP);
            end
        end

        [aggregation_min_abs_den_vector(t_Nr), ...
         aggregation_min_signed_den_vector(t_Nr), ...
         aggregation_negative_den_count_vector(t_Nr), ...
         aggregation_small_den_count_vector(t_Nr), ...
         aggregation_max_abs_phi_vector(t_Nr)] = ...
            consensus_aggregation_diagnostics(Xi_final, base_method, y_dim);

        if should_refresh_projection && mod(t_Nr-1,100) == 0
            if iscell(MaskedGP)
                masked_data_quantity = MaskedGP{1}.DataQuantity;
            else
                masked_data_quantity = MaskedGP.DataQuantity;
            end
            fprintf('Projection update at t = %.2f, MaskedGP data = %d.\n', ...
                t, masked_data_quantity);
        end
    end

    %% 9.3 Consensus law for tracking controller
    [phi_cell, ~, ~] = Manipulator_2D_2DoF_ConsensusLaw( ...
        vartheta_cell, x_tilde_cell, x_l_r, MultiAgentSystem, c, lambda_set, s_r_cell);

    %% 9.4 Prediction used by controller
    if is_ip_mode
        if isempty(MaskedGP)
            error('MaskedGP has not been initialized before control prediction.');
        end
        for AgentNr = 1:AgentQuantity
            if iscell(MaskedGP)
                current_masked_gp = MaskedGP{AgentNr};
            else
                current_masked_gp = MaskedGP;
            end
            if current_masked_gp.DataQuantity == 0
                error('MaskedGP has no inducing-point data before control prediction.');
            end
            [mu_hat, ~] = current_masked_gp.predict(x_all_matrix(:,AgentNr));
            mu_hat(~isfinite(mu_hat)) = 0;
            f_hat_matrix(:,AgentNr) = real(mu_hat);
        end

    elseif is_local_mode
        for AgentNr = 1:AgentQuantity
            [mu_hat, ~] = LocalGP_set{AgentNr}.predict(x_all_matrix(:,AgentNr));
            mu_hat(~isfinite(mu_hat)) = 0;
            f_hat_matrix(:,AgentNr) = real(mu_hat);
        end

    elseif is_without_gp_mode
        % Zero-mean / without-GP baseline: the plant still contains the true
        % unknown dynamics, but the controller applies no learned/exact
        % compensation for that unknown term.
        f_hat_matrix(:,:) = 0;

    elseif is_exact_mode
        for AgentNr = 1:AgentQuantity
            f_hat_matrix(:,AgentNr) = UnknownScale * ...
                Manipulator_2D_2DoF_UnknownDynamics(x_all_matrix(:,AgentNr));
        end
    end

    for AgentNr = 1:AgentQuantity
        f_true_matrix(:,AgentNr) = UnknownScale * ...
            Manipulator_2D_2DoF_UnknownDynamics(x_all_matrix(:,AgentNr));
    end
    prediction_error_norm_vector(t_Nr) = norm(f_true_matrix(:) - f_hat_matrix(:));
    f_hat_norm_vector(t_Nr) = norm(f_hat_matrix(:));
    f_hat_all_set(:,:,t_Nr) = f_hat_matrix;
    f_true_all_set(:,:,t_Nr) = f_true_matrix;

    if EnableProjectionDiagnostics && is_ip_mode
        for AgentNr = 1:AgentQuantity
            f_direct_matrix(:,AgentNr) = direct_aggregate_localgp_prediction( ...
                x_all_matrix(:,AgentNr), LocalGP_set, base_method, ...
                AgentQuantity, SigmaF);
            [mu_ip_cen, ~] = CenMaskedGP.predict(x_all_matrix(:,AgentNr));
            mu_ip_cen(~isfinite(mu_ip_cen)) = 0;
            f_ip_cen_matrix(:,AgentNr) = real(mu_ip_cen);
        end
        f_direct_all_set(:,:,t_Nr) = f_direct_matrix;
        f_ip_cen_all_set(:,:,t_Nr) = f_ip_cen_matrix;
        direct_prediction_error_norm_vector(t_Nr) = ...
            norm(f_true_matrix(:) - f_direct_matrix(:));
        ip_cen_prediction_error_norm_vector(t_Nr) = ...
            norm(f_true_matrix(:) - f_ip_cen_matrix(:));
        ip_projection_gap_norm_vector(t_Nr) = ...
            norm(f_hat_matrix(:) - f_direct_matrix(:));
    end

    %% 9.5 Control and physical dynamics
    u_cell = Manipulator_2D_2DoF_get_u_cell( ...
        x_all_cell, phi_cell, f_hat_matrix, L1, L2, m1, m2);

    [t_ode, x_all_temp] = ode45( ...
        @(current_time,current_state) Manipulator_2D_2DoF_MultiAgent_DynamicFunction( ...
            current_time, current_state, u_cell, L1, L2, m1, m2, UnknownScale), ...
        [t, t+t_step], x_all, opts);

    if isempty(t_ode) || t_ode(end) < t + t_step - 1e-10
        if ~isempty(SaveFolderName)
            debug_file = fullfile(SaveFolderName, [SaveFileName, '_ODE_FAIL_DEBUG.mat']);
            save(debug_file, 'CurrentMode', 't_Nr', 't', 'x_all', 'x_all_matrix', ...
                'u_cell', 'f_hat_matrix', 'f_true_matrix', 'P_inducing', 'Xi_final', ...
                'MaskedGP', 'LocalGP_set', 'InducingPoints_Coordinates', ...
                'prediction_error_norm_vector');
            error('ode45 failed at t = %.6f. Debug data saved to: %s', t, debug_file);
        else
            error('ode45 failed at t = %.6f.', t);
        end
    end

    x_all_next = x_all_temp(end,:)';
    x_all_set(:, t_Nr+1) = x_all_next;
    vartheta_all_set(:, t_Nr+1) = x_all_next - s_all_set(:,t_Nr+1) - ...
        kron(ones(AgentQuantity,1), xl_set(:,t_Nr+1));

    if mod(t_Nr,100) == 0
        elapsed_now = toc(elapsed_tic);
        progress_ratio = t_Nr / (T-1);
        estimated_total_time = elapsed_now / max(progress_ratio, eps);
        estimated_remaining_time = estimated_total_time - elapsed_now;
        %% Print Time
	    fprintf('t = %6.4f\n',t);

        fprintf(['t = %.2f / %.2f s, progress = %.1f%%, ', ...
                 'tracking error = %.4f, prediction error = %.4f, ', ...
                 'elapsed = %.1f min, remaining ≈ %.1f min\n'], ...
            t, t_end, 100*progress_ratio, ...
            TrackingError_vector(t_Nr), prediction_error_norm_vector(t_Nr), ...
            elapsed_now/60, estimated_remaining_time/60);
    end
end

TrackingError_vector(end) = norm(vartheta_all_set(:,end));
elapsed_time = toc(elapsed_tic);

% The two diagnostics required by the experiment goal.
AbsolutePChangeRMS = projection_absolute_change_rms_set;
FinitePChange = AbsolutePChangeRMS(isfinite(AbsolutePChangeRMS));
ComparisonCount = min(10,floor(numel(FinitePChange)/2));
if ComparisonCount > 0
    EarlyPChange = mean(FinitePChange(1:ComparisonCount));
    LatePChange = mean(FinitePChange(end-ComparisonCount+1:end));
else
    EarlyPChange = NaN;
    LatePChange = NaN;
end

CommunicationWindowCount = ceil(t_end);
DACBroadcastCountPerSecond = zeros(1,CommunicationWindowCount);
for WindowNr = 1:CommunicationWindowCount
    window_mask = t_set(1:T-1) >= WindowNr-1 & ...
        t_set(1:T-1) < WindowNr;
    DACBroadcastCountPerSecond(WindowNr) = ...
        sum(dac_broadcast_count_set(:,window_mask),'all')/AgentQuantity;
end

fprintf('\n==================================================\n');
fprintf('Mode: %s\n', CurrentMode);
fprintf('Total simulation time: %.2f s\n', elapsed_time);
fprintf('Final tracking error: %.6f\n', TrackingError_vector(end));
fprintf('Online updates: %.0f / agent\n', mean(online_update_count));
fprintf('UnknownScale: %.3f, DisturbanceScale: %.3f\n', UnknownScale, DisturbanceScale);
fprintf('Projection updates: %d\n', sum(projection_update_set));
if is_ip_mode
    fprintf('Consensus communication updates: %d\n', ...
        sum(consensus_communication_update_set));
end
if is_dac_mode
    fprintf('Mean absolute P-change RMS: early %.6g, late %.6g.\n', ...
        EarlyPChange,LatePChange);
    fprintf(['DAC packaged broadcasts per agent ', ...
        'in each one-second window:\n']);
    disp([(0:CommunicationWindowCount-1).', ...
        DACBroadcastCountPerSecond.']);
    fprintf(['DAC packaged broadcasts: %.0f total, ', ...
        '%.2f average per agent.\n'], ...
        sum(dac_broadcast_event_count_per_agent), ...
        mean(dac_broadcast_event_count_per_agent));
    fprintf(['DAC paper trigger (17): epsilon=%.4g. ', ...
        'No error-based terminal condition; %d DAC step per physical step.\n'], ...
        DACTriggerEpsilon,DACStepsPerPhysicalStep);
    completed_dac_mask = dac_iteration_count_set > 0;
    if any(completed_dac_mask)
        fprintf(['Last DAC step: %d fixed iteration, tracking error %.3e, ', ...
            'consensus disagreement %.3e.\n'], ...
            dac_iteration_count_set(find(completed_dac_mask,1,'last')), ...
            dac_tracking_error_set(find(completed_dac_mask,1,'last')), ...
            dac_consensus_disagreement_set(find(completed_dac_mask,1,'last')));
    end
elseif is_ac_mode
    fprintf('AC broadcasts: %.0f total, %.2f average per agent.\n', ...
        sum(ac_event_count_per_agent), ...
        mean(ac_event_count_per_agent));
    fprintf('AC online updates: %d. Communication steps: %d.\n', ...
        nnz(projection_update_set), nnz(ac_iteration_count_set > 0));
    fprintf(['AC periodic cached-broadcast trigger: sigma=%.4g; ', ...
        'Delta_c=%.4g s; kappa_P=%.4g; step=%.4g; ', ...
        '%d fixed rounds per 0.1 s P window; ', ...
        'no tolerance terminal.\n'], ...
        ACPeriodicSigma,ACCommunicationInterval,Kappa_P_AC, ...
        ACConsensusStep,ACFixedIterations);
end
if any(isfinite(prediction_error_norm_vector))
    [max_pred_err, idx] = max(prediction_error_norm_vector);
    fprintf('Max controller prediction error: %.6f at t=%.4f\n', max_pred_err, t_set(idx));
end
if any(isfinite(aggregation_min_abs_den_vector))
    [min_abs_den, idx] = min(aggregation_min_abs_den_vector, [], 'omitnan');
    fprintf('Min aggregation denominator magnitude: %.3e at t=%.4f\n', ...
        min_abs_den, t_set(idx));
    fprintf('Negative denominator samples: %.0f total; small denominator samples: %.0f total.\n', ...
        sum(aggregation_negative_den_count_vector, 'omitnan'), ...
        sum(aggregation_small_den_count_vector, 'omitnan'));
end
fprintf('==================================================\n');

if ~isempty(SaveFolderName) && ~isempty(SaveFileName)
    online_trigger_count = online_update_count;
    dac_broadcasts_per_agent = mean(dac_broadcast_event_count_per_agent);
    ac_broadcasts_per_agent = mean(ac_event_count_per_agent);
    ac_physical_triggers_per_agent = mean(sum(ac_physical_trigger_count_set,2));
    if ~exist('var_diag_min_set','var')
        % is_ip_mode was false, or the projection branch never ran --
        % keep the save() call below valid with empty placeholders.
        var_diag_min_set = [];
        var_diag_near_floor_set = [];
    end

    [~,~,save_ext] = fileparts(SaveFileName);
    if isempty(save_ext)
        save_file_name = [SaveFileName, '.mat'];
    else
        save_file_name = SaveFileName;
    end
    save(fullfile(SaveFolderName, save_file_name), ...
        'var_diag_min_set', 'var_diag_near_floor_set', ...
        't_set', 'TrackingError_vector', 'CurrentMode', 'use_formation', ...
        'x_all_set', 'vartheta_all_set', 'prediction_error_norm_vector', ...
        'f_hat_norm_vector', 'f_hat_all_set', 'f_true_all_set', ...
         'EnableProjectionDiagnostics', 'f_direct_all_set', ...
         'direct_prediction_error_norm_vector', ...
         'f_ip_cen_all_set', 'ip_cen_prediction_error_norm_vector', ...
         'ip_projection_gap_norm_vector', ...
        'aggregation_min_abs_den_vector', 'aggregation_min_signed_den_vector', ...
        'aggregation_negative_den_count_vector', 'aggregation_small_den_count_vector', ...
        'aggregation_max_abs_phi_vector', ...
        'online_update_count', 'online_trigger_count', ...
        'online_learning_trigger_count_set', 'online_learning_error_set', ...
        'projection_update_set', ...
        'consensus_communication_update_set', 'consensus_input_change_vector', ...
        'consensus_input_change_per_agent_set', ...
        'consensus_input_numerator_change_per_agent_set', ...
        'consensus_input_denominator_change_per_agent_set', ...
        'projection_step_change_per_agent_set', ...
        'projection_absolute_change_rms_set', 'AbsolutePChangeRMS', ...
        'projection_relative_change_per_agent_set', ...
        'projection_step_numerator_change_per_agent_set', ...
        'projection_step_denominator_change_per_agent_set', ...
        'consensus_input_trigger_count_set', ...
        'OfflineDataQuantity', 'is_offline_ip_mode', ...
        'IPOnlineUpdateInterval', 'IPOnlineUpdateStep', ...
        'OnlineAgentPolicy', ...
        'ACOnlineUpdateInterval', 'ACOnlineUpdateStep', ...
        'InitialSeedOverride', ...
         'UnknownScale', 'DisturbanceScale', 'SigmaL', ...
         'IPReconstructionSigmaN', ...
         'RBCMBetaMax', 'AggregationParameters', ...
         'InducingPoints_Coordinates', 'InducingPointSource', ...
         'InducingPointSeed', ...
         'P_inducing', 'Xi_final', 'elapsed_time', ...
        'NumInducingPoints', 'dac_broadcasts_per_agent', ...
        'ac_broadcasts_per_agent', 'dac_event_count_per_agent', ...
        'dac_broadcast_event_count_per_agent', ...
        'ac_event_count_per_agent', 'dac_broadcast_count_set', ...
        'dac_trigger_count_per_agent_point_set', ...
        'dac_inner_broadcast_count_set', ...
        'dac_iteration_count_set', 'dac_tracking_error_set', ...
        'dac_tracking_error_per_agent_set', ...
        'dac_consensus_disagreement_set', ...
        'EarlyPChange', 'LatePChange', ...
        'DACBroadcastCountPerSecond', ...
        'dac_inner_trigger_instance_set', ...
        'ac_broadcast_count_set', 'ac_inner_broadcast_count_set', ...
        'ac_physical_trigger_count_set', ...
        'ac_physical_triggers_per_agent', 'ac_iteration_count_set', ...
        'ac_consensus_disagreement_set', ...
        'ac_inner_trigger_instance_set', ...
        'ac_tracking_error_history_set', ...
        'EnableACTriggerDiagnostics', 'ac_trigger_measure_set', ...
        'ac_trigger_threshold_set', 'ACDetectorStoresSquared', ...
        'ConsensusInputTriggerTol', 'ACPeriodicSigma', ...
        'ACForceFullBroadcast', ...
         'DACTriggerEpsilon', 'Kappa_P_DAC', 'Kappa_P_AC', ...
         'ACConsensusStep', 'ACCommunicationInterval', ...
         'ACChecksPerPhysicalStep', 'ACTimingMode', ...
         'ACFixedIterations', ...
        'DACStepsPerPhysicalStep');
end

end

%% ========================================================================
% DAC equation (17): original agent-level packaged broadcast update
% ========================================================================
function [Xi_hat, triggered_points] = dac_equation17_broadcast_update( ...
    Xi, Xi_hat, L, epsilon_i, force_all)

[~,AgentQuantity,NumInducingPoints] = size(Xi);
triggered_points = false(AgentQuantity,NumInducingPoints);

for agent_i = 1:AgentQuantity
    neighbors = find(L(agent_i,:) < 0);
    degree_i = L(agent_i,agent_i);
    if degree_i <= 0 || isempty(neighbors)
        continue;
    end

    % One trigger broadcasts the agent's complete p-by-M matrix. Dividing
    % the Frobenius energies by M keeps epsilon comparable when M changes.
    local_error = Xi_hat(:,agent_i,:) - Xi(:,agent_i,:);
    local_error_sq = sum(local_error.^2,'all')/NumInducingPoints;

    neighbor_weights = reshape(-L(agent_i,neighbors),1,[],1);
    neighbor_gaps = Xi_hat(:,agent_i,:) - Xi_hat(:,neighbors,:);
    weighted_neighbor_gap_sq = sum( ...
        neighbor_weights.*neighbor_gaps.^2,'all')/NumInducingPoints;

    threshold = (weighted_neighbor_gap_sq+epsilon_i^2)/(4*degree_i);
    trigger_agent = force_all || local_error_sq > threshold;
    if trigger_agent
        triggered_points(agent_i,:) = true;
        Xi_hat(:,agent_i,:) = Xi(:,agent_i,:);
    end
end
end

% Rebuild one inducing-point GP for each agent from its current DAC output
% ========================================================================
function MaskedGP = build_agent_maskedgp_from_xi(Xi, method, AgentQuantity, M, p_dim, ...
    InducingPoints_Coordinates, SigmaF, SigmaL, x_dim, y_dim, MaskedGP)

reuse_factorization = iscell(MaskedGP) && numel(MaskedGP) == AgentQuantity && ...
    all(cellfun(@(gp) ~isempty(gp) && gp.DataQuantity == M, MaskedGP));

if ~reuse_factorization
    MaskedGP = cell(AgentQuantity,1);
end

for AgentNr = 1:AgentQuantity
    Xi_agent = Xi(:,AgentNr,:);
    phi_agent = decode_phi_from_xi(Xi_agent, method, AgentQuantity, M, p_dim, y_dim, SigmaF);

    if reuse_factorization
        MaskedGP{AgentNr} = update_maskedgp_targets(MaskedGP{AgentNr}, phi_agent);
    else
        reconstruction_sigma_n = get_env_positive_scalar( ...
            'CONTROL_IP_RECONSTRUCTION_SIGMA_N',1e-6);
        MaskedGP{AgentNr} = LocalGP_MultiOutput( ...
            x_dim,y_dim,M,reconstruction_sigma_n,SigmaF,SigmaL);
        MaskedGP{AgentNr}.add_Alldata(InducingPoints_Coordinates, phi_agent);
    end
end
end

%% ========================================================================
% Update only GP targets when inducing inputs and hyperparameters are fixed
% ========================================================================
function MaskedGP = update_maskedgp_targets(MaskedGP, phi)
M = MaskedGP.DataQuantity;
Y = phi.';
L = MaskedGP.L(1:M,1:M);
aux_alpha = L \ Y;
alpha = L' \ aux_alpha;

MaskedGP.Y(1:M,:) = Y;
MaskedGP.aux_alpha(1:M,:) = aux_alpha;
MaskedGP.alpha(1:M,:) = alpha;
end

% Dataset-style decoder: Xi_final -> phi -> shared MaskedGP
% ========================================================================
function MaskedGP = build_shared_maskedgp_from_xi(Xi_final, method, AgentQuantity, M, p_dim, ...
    InducingPoints_Coordinates, SigmaF, SigmaL, x_dim, y_dim)

phi = decode_phi_from_xi(Xi_final, method, AgentQuantity, M, p_dim, y_dim, SigmaF);
reconstruction_sigma_n = get_env_positive_scalar( ...
    'CONTROL_IP_RECONSTRUCTION_SIGMA_N',1e-6);
MaskedGP = LocalGP_MultiOutput( ...
    x_dim,y_dim,M,reconstruction_sigma_n,SigmaF,SigmaL);
MaskedGP.add_Alldata(InducingPoints_Coordinates, phi);
end

%% ========================================================================
% Exact centralized IP reference with reusable inducing-point factorization
% ========================================================================
function MaskedGP = build_or_update_shared_maskedgp_from_xi( ...
    Xi_final, method, AgentQuantity, M, p_dim, ...
    InducingPoints_Coordinates, SigmaF, SigmaL, x_dim, y_dim, MaskedGP)

phi = decode_phi_from_xi(Xi_final, method, AgentQuantity, M, p_dim, y_dim, SigmaF);
if isempty(MaskedGP) || MaskedGP.DataQuantity ~= M
    reconstruction_sigma_n = get_env_positive_scalar( ...
        'CONTROL_IP_RECONSTRUCTION_SIGMA_N',1e-6);
    MaskedGP = LocalGP_MultiOutput( ...
        x_dim,y_dim,M,reconstruction_sigma_n,SigmaF,SigmaL);
    MaskedGP.add_Alldata(InducingPoints_Coordinates, phi);
else
    MaskedGP = update_maskedgp_targets(MaskedGP, phi);
end
end

%% ========================================================================
% Decode one agent's DAC information vector into inducing-point targets
% ========================================================================
function phi = decode_phi_from_xi(Xi_final, method, AgentQuantity, M, p_dim, y_dim, SigmaF)
method = lower(method);
phi = zeros(y_dim, M);
prior_var = SigmaF^2;
aggregation_cfg = control_aggregation_parameters();
precision_floor = aggregation_cfg.precision_floor;
bcm_prior_scale = aggregation_cfg.bcm_prior_scale;

if p_dim < 2*y_dim
    error('p_dim=%d is too small for y_dim=%d.', p_dim, y_dim);
end

for d = 1:y_dim
    xi1 = squeeze(Xi_final(2*d-1, 1, :))';
    xi2 = squeeze(Xi_final(2*d,   1, :))';

    switch method
        case {'poe','gpoe','rbcm'}
            den = xi2;
            valid_den = den > precision_floor;
            phi(d,valid_den) = xi1(valid_den) ./ den(valid_den);

        case 'bcm'
            den_bcm = xi2;
            valid_bcm = den_bcm > precision_floor;
            phi(d,valid_bcm) = xi1(valid_bcm) ./ den_bcm(valid_bcm);

            % BCM can produce non-positive precision because of the
            % -(K-1)/prior_var correction.  When that happens during DAC
            % transients, fall back to the corresponding positive PoE
            % precision; BCM and PoE share the same numerator.
            den_poe = den_bcm + bcm_prior_scale * (AgentQuantity - 1) / prior_var;
            fallback = ~valid_bcm & den_poe > precision_floor;
            phi(d,fallback) = xi1(fallback) ./ den_poe(fallback);

        case 'moe'
            % Kept consistent with the dataset version.
            phi(d,:) = xi1 / AgentQuantity;

        otherwise
            error('Unknown aggregation method: %s.', method);
    end
end

phi(~isfinite(phi)) = 0;
phi_clip = 30;
phi = max(-phi_clip, min(phi_clip, phi));
end

%% ========================================================================
% Direct local-GP aggregation at the current query state for diagnostics
% ========================================================================
function f_direct = direct_aggregate_localgp_prediction( ...
    x_query, LocalGP_set, method, AgentQuantity, SigmaF)
method = lower(method);
y_dim = LocalGP_set{1}.y_dim;
prior_var = SigmaF^2;
aggregation_cfg = control_aggregation_parameters();
posterior_var_floor = aggregation_cfg.posterior_var_floor;
rbcm_beta_max = aggregation_cfg.rbcm_beta_max;
bcm_prior_scale = aggregation_cfg.bcm_prior_scale;
precision_floor = aggregation_cfg.precision_floor;

mu_set = nan(y_dim, AgentQuantity);
var_set = nan(y_dim, AgentQuantity);
for local_idx = 1:AgentQuantity
    [mu_i, var_i] = LocalGP_set{local_idx}.predict(x_query);
    mu_set(:,local_idx) = max(-30, min(30, mu_i));
    var_set(:,local_idx) = max(var_i, posterior_var_floor);
end

f_direct = zeros(y_dim,1);
for output_idx = 1:y_dim
    mu = mu_set(output_idx,:);
    var = var_set(output_idx,:);

    switch method
        case 'poe'
            num = sum(mu ./ var);
            den = sum(1 ./ var);

        case 'gpoe'
            beta = max(min(0.5*(log(prior_var)-log(var)), ...
                aggregation_cfg.gpoe_beta_max),eps);
            num = sum(beta .* mu ./ var);
            den = sum(beta ./ var);

        case 'moe'
            f_direct(output_idx) = mean(mu);
            continue;

        case 'bcm'
            num = sum(mu ./ var);
            den = sum(1 ./ var) - ...
                bcm_prior_scale * (AgentQuantity - 1) / prior_var;
            if den <= precision_floor
                den = sum(1 ./ var);
            end

        case 'rbcm'
            raw_beta = 0.5 * (log(prior_var) - log(var));
            beta = min(max(raw_beta, eps), rbcm_beta_max);
            num = sum(beta .* mu ./ var);
            den = sum(beta ./ var) + (1 - sum(beta)) / prior_var;

        otherwise
            error('Unknown aggregation method: %s', method);
    end

    if den > precision_floor
        f_direct(output_idx) = num / den;
    else
        f_direct(output_idx) = 0;
    end
end

f_direct(~isfinite(f_direct)) = 0;
f_direct = max(-30, min(30, f_direct));
end

function [min_abs_den, min_signed_den, negative_den_count, small_den_count, max_abs_phi] = ...
    consensus_aggregation_diagnostics(Xi_final, method, y_dim)
method = lower(method);

if isempty(Xi_final) || ~ismember(method, {'poe','gpoe','bcm','rbcm'})
    min_abs_den = NaN;
    min_signed_den = NaN;
    negative_den_count = NaN;
    small_den_count = NaN;
    max_abs_phi = NaN;
    return;
end

den_all = [];
phi_abs_all = [];
for d = 1:y_dim
    num_d = squeeze(Xi_final(2*d-1,:,:));
    den_d = squeeze(Xi_final(2*d,:,:));
    den_all = [den_all; den_d(:)]; %#ok<AGROW>

    valid_den = abs(den_d) > 1e-12;
    phi_d = nan(size(num_d));
    phi_d(valid_den) = num_d(valid_den) ./ den_d(valid_den);
    phi_abs_all = [phi_abs_all; abs(phi_d(:))]; %#ok<AGROW>
end

den_all = den_all(isfinite(den_all));
phi_abs_all = phi_abs_all(isfinite(phi_abs_all));

if isempty(den_all)
    min_abs_den = NaN;
    min_signed_den = NaN;
    negative_den_count = NaN;
    small_den_count = NaN;
else
    min_abs_den = min(abs(den_all));
    min_signed_den = min(den_all);
    negative_den_count = nnz(den_all < 0);
    small_den_count = nnz(abs(den_all) < 1e-4);
end

if isempty(phi_abs_all)
    max_abs_phi = NaN;
else
    max_abs_phi = max(phi_abs_all);
end
end

%% ========================================================================
% Apply graph Laplacian along the agent dimension of a 3D tensor
% ========================================================================
function L_X = laplacian_multiply_agent_dim_local(X, L)
% X has size p-by-N-by-M: information component, agent, inducing point.
% Put the N-agent dimension first, apply L to all M pages at once, and put
% the dimensions back. This is exactly sum_j L(i,j)*X_j for every component.
X_agent_first = permute(X,[2,1,3]);
L_X = permute(pagemtimes(L,X_agent_first),[2,1,3]);
end

function rms_values = rms_change_per_agent(delta_tensor, component_idx)
[~, agent_quantity, ~] = size(delta_tensor);
rms_values = nan(agent_quantity,1);
for agent_i = 1:agent_quantity
    delta_i = delta_tensor(component_idx,agent_i,:);
    rms_values(agent_i) = sqrt(mean(delta_i(:).^2, 'omitnan'));
end
end

function value = get_env_positive_scalar(name, default_value)
value = default_value;
raw = getenv(name);
if isempty(raw)
    return;
end
candidate = str2double(raw);
if isfinite(candidate) && candidate > 0
    value = candidate;
end
end

function value = get_env_positive_integer(name,default_value)
value = get_env_positive_scalar(name,default_value);
value = max(1,round(value));
end
