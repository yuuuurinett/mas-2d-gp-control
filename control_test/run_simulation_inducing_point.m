function [TrackingError_vector, t_set] = run_simulation_inducing_point(CurrentMode, SaveFolderName, SaveFileName, use_formation)
% Inducing-point online learning version requested by advisor:
%   1) Zero offline samples: each LocalGP starts from the prior.
%   2) No online data trigger: every agent adds one online sample at every time step.
%   3) Consensus communication ET is still kept:
%      - DAC modes: dynamic average consensus ET through gp_masked_aggregation_update.
%      - AC modes : static average consensus ET following the dataset IP-AC logic.
%
% Modes:
%   poe/gpoe/moe/bcm/rbcm       : IP-DAC + always-online data update
%   poe_ac/gpoe_ac/...          : IP-AC  + always-online data update
%   local                       : local GP + always-online data update
%   exact                       : exact unknown dynamics, no GP learning

if nargin < 4
    use_formation = true;
end

rng(0);

%% 1. System Parameters
SystemOrder = 2;
q_dim = 2;
x_dim = q_dim * SystemOrder;

m1 = 1;
m2 = 1;
L1 = 1;
L2 = 1;
g  = 9.8; 

AgentQuantity = 6;
LeaderQuantity = 1;

%% 2. Topology
MultiAgentSystem = Manipulator_2D_2DoF_SetMASTopology(AgentQuantity, LeaderQuantity);

L = MultiAgentSystem.Agent_Topology.LaplacianMatrix;
B = MultiAgentSystem.Agent_Leader_Topology.ConnectionMatrix(:,1);

%% 3. Controller Parameters
c = 10;
lambda_set = [1; 1];
Fl = 0.25;

lambda_n = lambda_set(end);
lambda_vector = lambda_set(1:SystemOrder-1);

Lambda = [zeros(SystemOrder-2,1), eye(SystemOrder-2); ...
          -lambda_set(1)/lambda_n, -lambda_set(2:end-1)/lambda_n];

% Pe & Qe
Qes = eye(SystemOrder-1);
Pes = care(Lambda, [], Qes);
Pe  = kron(Pes, eye(AgentQuantity));
Qe  = kron(Qes, eye(AgentQuantity));

% Pr & Qr
q  = (L + diag(B)) \ ones(AgentQuantity,1);
Pr = diag(1./q);
Qr = Pr*(L+diag(B)) + (L+diag(B))'*Pr;

% Pz & Qz
Pz = blkdiag(Pr, Pe);
eig_Pz = eig(Pz);
max_eig_Pz = max(eig_Pz);
min_eig_Pz = min(eig_Pz);

t_vec = [zeros(SystemOrder-2,1); 1/lambda_n];

Phi = Pr * kron(lambda_vector' * t_vec, eye(AgentQuantity));
Psi = Pr * kron(lambda_vector' * Lambda, eye(AgentQuantity)) + ...
      kron(t_vec' * Pes, eye(AgentQuantity));

Qz = [c*lambda_n*Qr - 2*Phi, -Psi; ...
      -Psi', Qe];

eig_Qz = eig(Qz);
min_eig_Qz = min(eig_Qz);

if all(real(eig_Qz) > 0) && all(real(eig(Lambda)) < 0)
    fprintf('The controller is stable!\n');
else
    error('Controller is not stable!');
end

xi = 2 * lambda_n / min_eig_Qz * norm(Pr * (L + diag(B)));
chi = sqrt((1 + norm([t_vec, Lambda])^2) * max_eig_Pz / min_eig_Pz) * ...
       norm(inv(L + diag(B)));

%% 4. Time
t_start = 0;
t_end = 10;
t_step = 0.01;

t_set = t_start:t_step:t_end;
T = numel(t_set);

%% 5. Reference Trajectory
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

%% 6. Local GPs
SigmaF = 1;
SigmaL = 0.5 * ones(x_dim,1);

GP_tau = 1e-8;
GP_delta = 0.1;
y_dim = q_dim;
LocalGP_Quantity = AgentQuantity;

DomainScale = 1.5;
X_min = DomainScale * [-1,-1,-1,-1];
X_max = DomainScale * [ 1, 1, 1, 1];

MaxDataQuantity_set = 600 * ones(AgentQuantity,1);
SigmaN_set = 0.01 * ones(AgentQuantity,1);

dac_methods = {'poe','gpoe','moe','bcm','rbcm'};
ac_methods  = {'poe_ac','gpoe_ac','moe_ac','bcm_ac','rbcm_ac'};

mode_lower = lower(CurrentMode);

% Current advisor-requested online-learning setting:
% start from zero offline samples, then add online data during the closed-loop simulation.
% This makes the comparison focus on online learning instead of offline pretraining.
OfflineDataQuantity_set = zeros(AgentQuantity, 1);

LocalGP_set = cell(LocalGP_Quantity, 1);

for LocalGP_Nr = 1:LocalGP_Quantity
    MaxDataQuantity = MaxDataQuantity_set(LocalGP_Nr);
    OfflineDataQuantity = OfflineDataQuantity_set(LocalGP_Nr);
    SigmaN = SigmaN_set(LocalGP_Nr);

    LocalGP_set{LocalGP_Nr} = LocalGP_MultiOutput( ...
        x_dim, y_dim, MaxDataQuantity, SigmaN, SigmaF, SigmaL);

    if OfflineDataQuantity > 0
        X_in = 2*(rand(x_dim, OfflineDataQuantity)-0.5)*DomainScale;
        Y_in = Manipulator_2D_2DoF_UnknownDynamics(X_in);
        Y_in = Y_in + SigmaN * randn(size(Y_in));

        LocalGP_set{LocalGP_Nr}.add_Alldata(X_in, Y_in);
    end
    LocalGP_set{LocalGP_Nr}.tau   = GP_tau;
    LocalGP_set{LocalGP_Nr}.delta = GP_delta;
    LocalGP_set{LocalGP_Nr}.xMax  = X_max;
    LocalGP_set{LocalGP_Nr}.xMin  = X_min;
end

%% 7. Online Data Update Setting
% Advisor-requested setting: no online trigger.
% Every GP-based mode adds the current sample at every time step.
% Therefore eta-bound computation and DistributedET are not used in this version.
online_learning_modes = [dac_methods, ac_methods, {'local'}];

%% 9. Inducing Points & Aggregation Init
Kappa_P = 1;
NumInducingPoints = 400;
InducingPoints_Coordinates = 2*DomainScale*rand(x_dim, NumInducingPoints) - DomainScale;
L_lap = MultiAgentSystem.Agent_Topology.LaplacianMatrix;

% Dataset-style consensus ET parameters.
N_degree = sum(L_lap < 0, 2);
N_max = max(N_degree);
a_param = 0.5 / N_max;
sigma_i_dac = 0.5;
sigma_i_ac  = 0.5;

MaskedGP = [];
P_inducing = [];
p_dim = [];

% DAC state variables.
Zeta_vector_inducing = [];
Zeta_last_trigger = [];
dac_total_trigger_count = zeros(AgentQuantity, 1);
neighbor_count_per_agent = sum(L_lap < 0, 2);

% AC state variables.
Xi_ac = [];
Xi_last_trigger_ac = [];

% AC communication statistics are separated into two levels:
%   1) broadcast count: one communication event per agent and time step
%      if at least one inducing point triggers;
%   2) point count: point-wise trigger events for fair comparison with DAC.
ac_total_broadcast_count = zeros(AgentQuantity, 1);
ac_total_point_count     = zeros(AgentQuantity, NumInducingPoints);

% Backward-compatible alias used by some batch scripts.
ac_total_trigger_count = ac_total_broadcast_count;

if ismember(mode_lower, dac_methods)
    base_method = mode_lower;

    warning_state_gp_init = warning('off','all');
    [P_inducing, p_dim] = gp_masked_aggregation_init( ...
        LocalGP_set, AgentQuantity, NumInducingPoints, ...
        InducingPoints_Coordinates, base_method);
    warning(warning_state_gp_init);

    Zeta_vector_inducing = zeros(p_dim, AgentQuantity, NumInducingPoints);
    Zeta_last_trigger = P_inducing;

    [MaskedGP, Zeta_vector_inducing, Zeta_last_trigger, ~] = ...
        gp_masked_aggregation_update( ...
        P_inducing, Zeta_vector_inducing, L_lap, Kappa_P, AgentQuantity, ...
        NumInducingPoints, 0, InducingPoints_Coordinates, ...
        SigmaF, SigmaL, x_dim, base_method, p_dim, ...
        Zeta_last_trigger, neighbor_count_per_agent);

elseif ismember(mode_lower, ac_methods)
    base_method = strrep(mode_lower, '_ac', '');

    warning_state_gp_init = warning('off','all');
    [P_inducing, p_dim] = gp_masked_aggregation_init( ...
        LocalGP_set, AgentQuantity, NumInducingPoints, ...
        InducingPoints_Coordinates, base_method);
    warning(warning_state_gp_init);

    Xi_ac = P_inducing;
    Xi_last_trigger_ac = P_inducing;

    [MaskedGP, Xi_ac, Xi_last_trigger_ac, ~, ~] = ...
        gp_masked_aggregation_ac_et_update( ...
        Xi_ac, Xi_last_trigger_ac, L_lap, Kappa_P, AgentQuantity, ...
        NumInducingPoints, 0, InducingPoints_Coordinates, ...
        SigmaF, SigmaL, x_dim, base_method, p_dim, ...
        N_degree, a_param, sigma_i_ac);

elseif strcmpi(mode_lower, 'local') || strcmpi(mode_lower, 'offline') || strcmpi(mode_lower, 'exact')
    base_method = mode_lower;
else
    error('Unknown CurrentMode: %s', CurrentMode);
end

%% 10. Initial State
rng(42);

x_all = rand(x_dim*AgentQuantity, 1);

x_all_set = nan(x_dim*AgentQuantity, T);
x_all_set(:,1) = x_all;

vartheta_all_set = nan(x_dim*AgentQuantity, T);
vartheta_all_set(:,1) = x_all - s_all_set(:,1) - ...
    kron(ones(AgentQuantity,1), xl_set(:,1));

f_hat_matrix  = zeros(y_dim, AgentQuantity);
f_true_matrix = zeros(y_dim, AgentQuantity);

TrackingError_vector = zeros(1, T);

f_hat_all_set  = nan(y_dim, AgentQuantity, T);
f_true_all_set = nan(y_dim, AgentQuantity, T);

% Online data update logging variables.
% Since online data update is performed at every step, these variables are
% kept only as sanity-check logs and are not printed as trigger statistics.
online_update_set   = zeros(AgentQuantity, T);
online_update_count = zeros(AgentQuantity, 1);

%% 11. Control Loop
opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
tic;

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

    %% 11.1 Consensus law for controller
    [phi_cell, r_matrix, e_cell] = Manipulator_2D_2DoF_ConsensusLaw( ...
        vartheta_cell, x_tilde_cell, x_l_r, MultiAgentSystem, c, lambda_set, s_r_cell);

    %% 11.2 Prediction using current model
    if ismember(mode_lower, dac_methods) || ismember(mode_lower, ac_methods)
        for n = 1:AgentQuantity
            [mu_hat,~] = MaskedGP{n}.predict(x_all_matrix(:,n));
            f_hat_matrix(:,n) = max(-30, min(30, mu_hat));
        end

    elseif strcmpi(mode_lower, 'local') || strcmpi(mode_lower, 'offline')
        for n = 1:AgentQuantity
            [mu_hat,~] = LocalGP_set{n}.predict(x_all_matrix(:,n));
            f_hat_matrix(:,n) = max(-30, min(30, mu_hat));
        end

    elseif strcmpi(mode_lower, 'exact')
        for n = 1:AgentQuantity
            f_hat_matrix(:,n) = Manipulator_2D_2DoF_UnknownDynamics(x_all_matrix(:,n));
        end
    end

    %% 11.3 Record prediction and true dynamics
    for n = 1:AgentQuantity
        f_true_matrix(:,n) = Manipulator_2D_2DoF_UnknownDynamics(x_all_matrix(:,n));
    end

    f_hat_all_set(:,:,t_Nr)  = f_hat_matrix;
    f_true_all_set(:,:,t_Nr) = f_true_matrix;

    %% 11.4 Control input and system simulation
    u_cell = Manipulator_2D_2DoF_get_u_cell( ...
        x_all_cell, phi_cell, f_hat_matrix, L1, L2, m1, m2);

    [~, x_all_temp] = ode45( ...
        @(t,x) Manipulator_2D_2DoF_MultiAgent_DynamicFunction(t, x, u_cell, L1, L2, m1, m2), ...
        [t, t+t_step], x_all, opts);

    x_all_next = x_all_temp(end,:)';
    x_all_set(:, t_Nr+1) = x_all_next;

    vartheta_all_set(:, t_Nr+1) = x_all_next - s_all_set(:,t_Nr+1) - ...
        kron(ones(AgentQuantity,1), xl_set(:,t_Nr+1));

    %% 11.5 Online data update layer: always add one sample per agent
    % No online event trigger is used in this version.
    % Every GP-based agent adds the current sample (x_i, y_i) to its LocalGP at every step.
    any_online_update = false;
    updated_agents = [];

    if ismember(mode_lower, online_learning_modes)
        online_update_set(:, t_Nr) = 1;   % kept as a data-update indicator for logging

        for AgentNr = 1:AgentQuantity
            x_i = x_all_matrix(:, AgentNr);

            y_i = Manipulator_2D_2DoF_UnknownDynamics(x_i) + ...
                  LocalGP_set{AgentNr}.SigmaN * randn(y_dim, 1);

            if LocalGP_set{AgentNr}.DataQuantity >= LocalGP_set{AgentNr}.MaxDataQuantity
                LocalGP_set{AgentNr}.downdateParam(1);
            end

            LocalGP_set{AgentNr}.addPoint(x_i, y_i);

            online_update_count(AgentNr) = online_update_count(AgentNr) + 1;
            any_online_update = true;
            updated_agents = [updated_agents; AgentNr]; %#ok<AGROW>
        end
    end

    %% 11.6 Aggregation update layer for next step
    if ismember(mode_lower, dac_methods)
        % If online learning updated LocalGP, refresh only changed agents in P_inducing.
        if any_online_update
            for updatedAgentIdx = 1:numel(updated_agents)
                UpdatedAgentNr = updated_agents(updatedAgentIdx);

                P_inducing = recompute_P_single_agent( ...
                    P_inducing, UpdatedAgentNr, LocalGP_set, ...
                    NumInducingPoints, InducingPoints_Coordinates, ...
                    AgentQuantity, base_method);
            end
        end

        % DAC consensus ET: dynamic average consensus layer.
        [MaskedGP_new, Zeta_vector_inducing, Zeta_last_trigger, dac_step_triggers] = ...
            gp_masked_aggregation_update( ...
            P_inducing, Zeta_vector_inducing, L_lap, Kappa_P, AgentQuantity, ...
            NumInducingPoints, t_step, InducingPoints_Coordinates, ...
            SigmaF, SigmaL, x_dim, base_method, p_dim, ...
            Zeta_last_trigger, neighbor_count_per_agent);

        if any(dac_step_triggers)
            MaskedGP = MaskedGP_new;
        end

        dac_total_trigger_count = dac_total_trigger_count + dac_step_triggers;

    elseif ismember(mode_lower, ac_methods)
        % AC consensus ET: static average consensus layer, dataset-style.
        % Static AC assumes a fixed input Pi. If online data changes Pi,
        % refresh the current AC state Xi_ac with the new P_inducing.
        % The last transmitted state Xi_last_trigger_ac is NOT refreshed for free.
        if any_online_update
            for updatedAgentIdx = 1:numel(updated_agents)
                UpdatedAgentNr = updated_agents(updatedAgentIdx);

                P_inducing = recompute_P_single_agent( ...
                    P_inducing, UpdatedAgentNr, LocalGP_set, ...
                    NumInducingPoints, InducingPoints_Coordinates, ...
                    AgentQuantity, base_method);
            end

            Xi_ac = P_inducing;
            % Do NOT reset Xi_last_trigger_ac here.
            % Xi_last_trigger_ac stores the last actually transmitted AC state
            % and should only be updated when the AC communication trigger fires.
        end

        [MaskedGP, Xi_ac, Xi_last_trigger_ac, ...
         ac_step_broadcast_triggers, ac_step_point_triggers] = ...
            gp_masked_aggregation_ac_et_update( ...
            Xi_ac, Xi_last_trigger_ac, L_lap, Kappa_P, AgentQuantity, ...
            NumInducingPoints, t_step, InducingPoints_Coordinates, ...
            SigmaF, SigmaL, x_dim, base_method, p_dim, ...
            N_degree, a_param, sigma_i_ac);

        ac_total_broadcast_count = ...
            ac_total_broadcast_count + ac_step_broadcast_triggers;

        ac_total_point_count = ...
            ac_total_point_count + ac_step_point_triggers;

        % Backward-compatible alias: old scripts interpret this as AC / agent.
        ac_total_trigger_count = ac_total_broadcast_count;
    end

    % fprintf('t = %6.4f\n', t);
end

TrackingError_vector(end) = norm(vartheta_all_set(:,end));

%% 12. Final Record and Statistics

TrackingError_vector(end) = norm(vartheta_all_set(:,end));

x_all_matrix_end = reshape(x_all_set(:,end), x_dim, AgentQuantity);

for n = 1:AgentQuantity
    f_true_all_set(:,n,end) = Manipulator_2D_2DoF_UnknownDynamics(x_all_matrix_end(:,n));
    f_hat_all_set(:,n,end)  = f_hat_matrix(:,n);
end

elapsed_time = toc;

fprintf('\n==================================================\n');
fprintf('Mode: %s\n', CurrentMode);
fprintf('Formation: %d\n', use_formation);
fprintf('Total simulation time: %.2f s\n', elapsed_time);
fprintf('Final tracking error: %.6f\n', TrackingError_vector(end));

dac_trigger_count_per_agent_point = 0;
ac_trigger_count_per_agent        = 0;
ac_trigger_count_per_agent_point  = 0;

if ismember(lower(CurrentMode), dac_methods)

    dac_trigger_count_per_agent_point = ...
        mean(dac_total_trigger_count) / NumInducingPoints;

    fprintf('\nDAC consensus ET:\n');
    fprintf('  Average triggers: %.4f / agent / inducing point\n', ...
        dac_trigger_count_per_agent_point);

elseif ismember(lower(CurrentMode), ac_methods)

    ac_trigger_count_per_agent = mean(ac_total_broadcast_count);

    ac_trigger_count_per_agent_point = ...
        mean(sum(ac_total_point_count, 2)) / NumInducingPoints;

    fprintf('\nAC consensus ET:\n');
    fprintf('  Broadcast triggers: %.2f / agent\n', ...
        ac_trigger_count_per_agent);
    fprintf('  Point triggers: %.4f / agent / inducing point\n', ...
        ac_trigger_count_per_agent_point);
end

fprintf('==================================================\n');
%% 13. Save
if nargin >= 3
    if ~exist(SaveFolderName,'dir')
        mkdir(SaveFolderName);
    end

    save(fullfile(SaveFolderName,[SaveFileName,'.mat']), ...
        't_set', ...
        'TrackingError_vector', ...
        'CurrentMode', ...
        'use_formation', ...
        'f_hat_all_set', ...
        'f_true_all_set', ...
        'vartheta_all_set', ...
        'online_update_set', ...
        'online_update_count', ...
        'dac_total_trigger_count', ...
        'dac_trigger_count_per_agent_point', ...
        'ac_total_trigger_count', ...
        'ac_total_broadcast_count', ...
        'ac_total_point_count', ...
        'ac_trigger_count_per_agent', ...
        'ac_trigger_count_per_agent_point', ...
        'NumInducingPoints', ...
        'Kappa_P', ...
        'elapsed_time', ...
        'OfflineDataQuantity_set');
end

end

%% ========================================================================
%  Local helper: dataset-style IP-AC consensus ET update
%  ========================================================================
function [MaskedGP, Xi, Xi_last_trigger, broadcast_trigger_count, point_trigger_count] = gp_masked_aggregation_ac_et_update( ...
    Xi, Xi_last_trigger, L, Kappa_P, AgentQuantity, NumInducingPoints, TimeStep, ...
    InducingPoints_Coordinates, SigmaF, SigmaL, x_dim, method, p_dim, ...
    N_degree, a_param, sigma_i_ac)

M = NumInducingPoints;
y_dim = 2;
broadcast_trigger_count = zeros(AgentQuantity, 1);
point_trigger_count     = zeros(AgentQuantity, M);

if TimeStep > 0
    L_Xi_hat = laplacian_multiply_agent_dim_local(Xi_last_trigger, L);
    Xi = Xi - TimeStep * Kappa_P * L_Xi_hat;

    for agent_i = 1:AgentQuantity
        neighbor_i = (L(agent_i,:) < 0);
        N_i = N_degree(agent_i);

        if N_i <= 0
            continue;
        end

        coeff_i = sigma_i_ac * a_param * (1 - a_param*N_i) / N_i;

        E_i = Xi_last_trigger(:,agent_i,:) - Xi(:,agent_i,:);
        e_norm_sq = squeeze(sum(E_i.^2, 1));

        Z_i = N_i*Xi(:,agent_i,:) - sum(Xi(:,neighbor_i,:), 2);
        z_norm_sq = squeeze(sum(Z_i.^2, 1));

        trigger_idx = (e_norm_sq(:).' > coeff_i * z_norm_sq(:).');

        if any(trigger_idx)
            Xi_last_trigger(:,agent_i,trigger_idx) = Xi(:,agent_i,trigger_idx);

            % One broadcast event for this agent at this time step.
            broadcast_trigger_count(agent_i) = 1;

            % Point-wise record, useful for DAC-like normalized statistics.
            point_trigger_count(agent_i,trigger_idx) = 1;
        end
    end
end

MaskedGP = build_maskedgp_from_consensus_state( ...
    Xi, method, AgentQuantity, M, p_dim, InducingPoints_Coordinates, ...
    SigmaF, SigmaL, x_dim, y_dim);

end

%% ========================================================================
%  Local helper: build MaskedGP cell from consensus state Xi_final
%  ========================================================================
function MaskedGP = build_maskedgp_from_consensus_state( ...
    Xi_final, method, AgentQuantity, M, p_dim, InducingPoints_Coordinates, ...
    SigmaF, SigmaL, x_dim, y_dim)

prior_var = SigmaF^2; %#ok<NASGU>
method = lower(method);

num1 = squeeze(Xi_final(1, :, :));
aux1 = squeeze(Xi_final(2, :, :));
num2 = squeeze(Xi_final(3, :, :));
aux2 = squeeze(Xi_final(4, :, :));

switch method
    case {'poe','gpoe','bcm','rbcm'}
        phi1 = zeros(size(num1));
        phi2 = zeros(size(num2));
        mask1 = abs(aux1) > 1e-8;
        mask2 = abs(aux2) > 1e-8;
        phi1(mask1) = num1(mask1) ./ aux1(mask1);
        phi2(mask2) = num2(mask2) ./ aux2(mask2);

    case 'moe'
        phi1 = num1 / AgentQuantity;
        phi2 = num2 / AgentQuantity;

    otherwise
        error('Unknown aggregation method: %s', method);
end

phi1(~isfinite(phi1)) = 0;
phi2(~isfinite(phi2)) = 0;

MaskedGP = cell(AgentQuantity, 1);
for AgentNr = 1:AgentQuantity
    Y_agent = [phi1(AgentNr, :); phi2(AgentNr, :)];
    MaskedGP{AgentNr} = LocalGP_MultiOutput(x_dim, y_dim, M, 1e-4, SigmaF, SigmaL);
    MaskedGP{AgentNr}.add_Alldata(InducingPoints_Coordinates, Y_agent);
end

end

%% ========================================================================
%  Local helper: apply graph Laplacian along agent dimension
%  ========================================================================
function L_X = laplacian_multiply_agent_dim_local(X, L)
[p_dim, agent_quantity, num_points] = size(X);
X_agent_first = permute(X, [2, 1, 3]);
X_flat = reshape(X_agent_first, agent_quantity, []);
L_X_flat = L * X_flat;
L_X_agent_first = reshape(L_X_flat, agent_quantity, p_dim, num_points);
L_X = permute(L_X_agent_first, [2, 1, 3]);
end
