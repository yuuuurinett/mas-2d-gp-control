function [TrackingError_vector, t_set] = run_simulation_inducing_point(CurrentMode, SaveFolderName, SaveFileName, use_formation)
% Pure inducing-point online-learning simulation.
%
% Purpose of this version:
%   Keep only the core algorithmic logic for advisor discussion.
%
% Core logic:
%   1) Offline training data = 0.
%   2) At every sampling instant t_k, each agent adds one online sample
%          (x_i(t_k), y_i(t_k)),  y_i = UnknownScale*UnknownDynamics(x_i) + scaled noise.
%      There is no lower-level learning event trigger.
%   3) Every ProjectionUpdatePeriod steps, all LocalGPs are projected onto
%      a fixed inducing-point set.
%   4) The inducing-point statistics are aggregated by either IP-DAC or
%      IP-AC using Kia et al. style asynchronous distributed event-triggered
%      communication for connected undirected graphs.
%   5) A shared MaskedGP is rebuilt from the final consensus state.
%   6) For inducing-point modes, the controller always uses this shared
%      inducing-point MaskedGP prediction. There is no shadow/local switch.
%
% Supported modes:
%   poe, gpoe, moe, bcm, rbcm                 : IP-DAC
%   poe_ac, gpoe_ac, moe_ac, bcm_ac, rbcm_ac : IP-AC
%   local                                    : LocalGP only baseline
%   exact                                    : exact unknown dynamics baseline
%
% Note:
%   This file assumes that gp_masked_aggregation_init and the decoder below
%   use the same row layout for P / Xi_final as in the dataset code.

fprintf('Using PURE_SCALED_UNKNOWN_KIA_IPCONTROL.\n');

if nargin < 4
    use_formation = true;
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
ProjectionUpdatePeriod = 50;     % 50 steps with dt=0.01 -> refresh every 0.5 s
ConsensusMaxIter = 3000;
ConsensusTol = 1e-5;

% Scaling test suggested for numerical robustness.
% UnknownScale scales the unknown dynamics used for GP labels and prediction-error evaluation
% in this script. If your physical dynamics function internally calls
% Manipulator_2D_2DoF_UnknownDynamics, apply the same scale inside that function as well
% for a fully consistent scaled-plant experiment.
UnknownScale = 0.1;
DisturbanceScale = 0.1;

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
N_degree = sum(L_lap < 0, 2);

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
% For advisor discussion, keep this normal 10 s horizon. For debugging, you
% can temporarily change it to 1 or 2.
t_start = 0;
t_end = 10;
t_step = 0.01;
t_set = t_start:t_step:t_end;
T = numel(t_set);

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

%% 6. LocalGPs: zero offline data
SigmaF = 1;
SigmaL = 0.5 * ones(x_dim,1);
SigmaN = 0.01;
GP_tau = 1e-8;
GP_delta = 0.1;
DomainScale = 1.5;
X_min = DomainScale * [-1,-1,-1,-1];
X_max = DomainScale * [ 1, 1, 1, 1];

MaxDataQuantity = 600;
OfflineDataQuantity = 0;

LocalGP_set = cell(AgentQuantity, 1);
for AgentNr = 1:AgentQuantity
    LocalGP_set{AgentNr} = LocalGP_MultiOutput(x_dim, y_dim, MaxDataQuantity, ...
        SigmaN, SigmaF, SigmaL);
    LocalGP_set{AgentNr}.tau   = GP_tau;
    LocalGP_set{AgentNr}.delta = GP_delta;
    LocalGP_set{AgentNr}.xMax  = X_max;
    LocalGP_set{AgentNr}.xMin  = X_min;
end

fprintf('Offline data: %d per agent. Online update: every sampling step.\n', OfflineDataQuantity);
fprintf('UnknownScale = %.3f, DisturbanceScale = %.3f.\n', UnknownScale, DisturbanceScale);

%% 7. Mode and inducing-point setup
mode_lower = lower(CurrentMode);
dac_methods = {'poe','gpoe','moe','bcm','rbcm'};
ac_methods  = {'poe_ac','gpoe_ac','moe_ac','bcm_ac','rbcm_ac'};

is_dac_mode = ismember(mode_lower, dac_methods);
is_ac_mode  = ismember(mode_lower, ac_methods);
is_ip_mode  = is_dac_mode || is_ac_mode;
is_local_mode = strcmpi(mode_lower, 'local');
is_exact_mode = strcmpi(mode_lower, 'exact');

if is_dac_mode
    base_method = mode_lower;
elseif is_ac_mode
    base_method = strrep(mode_lower, '_ac', '');
elseif is_local_mode || is_exact_mode
    base_method = mode_lower;
else
    error('Unknown mode: %s', CurrentMode);
end

Kappa_P = 1;
NumInducingPoints = 400;
InducingPoints_Coordinates = 2*DomainScale*rand(x_dim, NumInducingPoints) - DomainScale;

% Kia et al. connected-undirected-graph trigger parameter.
epsilon_i = 0.04 * ones(AgentQuantity, 1);

%% 8. Initial state
rng(42);
x_all = rand(x_dim*AgentQuantity, 1);
x_all_set = nan(x_dim*AgentQuantity, T);
x_all_set(:,1) = x_all;

vartheta_all_set = nan(x_dim*AgentQuantity, T);
vartheta_all_set(:,1) = x_all - s_all_set(:,1) - kron(ones(AgentQuantity,1), xl_set(:,1));

TrackingError_vector = zeros(1,T);
f_hat_matrix = zeros(y_dim, AgentQuantity);
f_true_matrix = zeros(y_dim, AgentQuantity);

prediction_error_norm_vector = nan(1,T);
online_update_count = zeros(AgentQuantity,1);
projection_update_set = zeros(1,T);

MaskedGP = [];
Xi_final = [];
P_inducing = [];
p_dim = [];

fprintf('Mode: %s. Projection refresh every %d steps.\n', CurrentMode, ProjectionUpdatePeriod);
if is_ip_mode
    fprintf('Control source: shared inducing-point MaskedGP.\n');
elseif is_local_mode
    fprintf('Control source: LocalGP baseline.\n');
elseif is_exact_mode
    fprintf('Control source: exact unknown dynamics.\n');
end

%% 9. Main loop
opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
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

    %% 9.1 Time-triggered LocalGP online update
    if is_ip_mode || is_local_mode
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

    %% 9.2 Inducing-point projection and Kia ET consensus
    if is_ip_mode && (t_Nr == 1 || mod(t_Nr-1, ProjectionUpdatePeriod) == 0)
        projection_update_set(t_Nr) = 1;

        [P_inducing, p_dim] = gp_masked_aggregation_init( ...
            LocalGP_set, AgentQuantity, NumInducingPoints, ...
            InducingPoints_Coordinates, base_method);

        if is_dac_mode
            Xi_final = ip_dac_consensus_kia(P_inducing, L_lap, Kappa_P, t_step, ...
                ConsensusMaxIter, ConsensusTol, N_degree, epsilon_i);
        else
            Xi_final = ip_ac_consensus_kia(P_inducing, L_lap, Kappa_P, t_step, ...
                ConsensusMaxIter, ConsensusTol, N_degree, epsilon_i);
        end

        MaskedGP = build_shared_maskedgp_from_xi(Xi_final, base_method, AgentQuantity, ...
            NumInducingPoints, p_dim, InducingPoints_Coordinates, SigmaF, SigmaL, x_dim, y_dim);

        fprintf('Projection update at t = %.2f, MaskedGP data = %d.\n', t, MaskedGP.DataQuantity);
    end

    %% 9.3 Consensus law for tracking controller
    [phi_cell, ~, ~] = Manipulator_2D_2DoF_ConsensusLaw( ...
        vartheta_cell, x_tilde_cell, x_l_r, MultiAgentSystem, c, lambda_set, s_r_cell);

    %% 9.4 Prediction used by controller
    if is_ip_mode
        if isempty(MaskedGP) || MaskedGP.DataQuantity == 0
            error('MaskedGP has not been initialized before control prediction.');
        end
        for AgentNr = 1:AgentQuantity
            [mu_hat, ~] = MaskedGP.predict(x_all_matrix(:,AgentNr));
            mu_hat(~isfinite(mu_hat)) = 0;
            f_hat_matrix(:,AgentNr) = real(mu_hat);
        end

    elseif is_local_mode
        for AgentNr = 1:AgentQuantity
            [mu_hat, ~] = LocalGP_set{AgentNr}.predict(x_all_matrix(:,AgentNr));
            mu_hat(~isfinite(mu_hat)) = 0;
            f_hat_matrix(:,AgentNr) = real(mu_hat);
        end

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

    %% 9.5 Control and physical dynamics
    u_cell = Manipulator_2D_2DoF_get_u_cell( ...
        x_all_cell, phi_cell, f_hat_matrix, L1, L2, m1, m2);

    [t_ode, x_all_temp] = ode45( ...
        @(current_time,current_state) Manipulator_2D_2DoF_MultiAgent_DynamicFunction( ...
            current_time, current_state, u_cell, L1, L2, m1, m2), ...
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

fprintf('\n==================================================\n');
fprintf('Mode: %s\n', CurrentMode);
fprintf('Total simulation time: %.2f s\n', elapsed_time);
fprintf('Final tracking error: %.6f\n', TrackingError_vector(end));
fprintf('Online updates: %.0f / agent\n', mean(online_update_count));
fprintf('UnknownScale: %.3f, DisturbanceScale: %.3f\n', UnknownScale, DisturbanceScale);
fprintf('Projection updates: %d\n', sum(projection_update_set));
if any(isfinite(prediction_error_norm_vector))
    [max_pred_err, idx] = max(prediction_error_norm_vector);
    fprintf('Max controller prediction error: %.6f at t=%.4f\n', max_pred_err, t_set(idx));
end
fprintf('==================================================\n');

if ~isempty(SaveFolderName) && ~isempty(SaveFileName)
    save(fullfile(SaveFolderName, [SaveFileName, '.mat']), ...
        't_set', 'TrackingError_vector', 'CurrentMode', 'use_formation', ...
        'x_all_set', 'vartheta_all_set', 'prediction_error_norm_vector', ...
        'online_update_count', 'projection_update_set', ...
        'UnknownScale', 'DisturbanceScale', ...
        'InducingPoints_Coordinates', 'P_inducing', 'Xi_final', 'elapsed_time');
end

end

%% ========================================================================
% IP-DAC with Kia-style event-triggered communication
% ========================================================================
function Xi_final = ip_dac_consensus_kia(Pi, L, Kappa_P, t_step, max_iter, tol, N_degree, epsilon_i)
[p_dim, AgentQuantity, NumInducingPoints] = size(Pi);
Zeta = zeros(p_dim, AgentQuantity, NumInducingPoints);
Xi_hat = Pi;

for iter = 1:max_iter
    Zeta_prev = Zeta;

    L_Xi_hat = laplacian_multiply_agent_dim_local(Xi_hat, L);
    Zeta = Zeta + t_step * Kappa_P * L_Xi_hat;
    Xi = Pi - Zeta;

    Xi_hat = kia_trigger_update(Xi, Xi_hat, L, N_degree, epsilon_i, iter);

    if max(abs(Zeta(:) - Zeta_prev(:))) < tol
        break;
    end
end

Xi_final = Pi - Zeta;
end

%% ========================================================================
% IP-AC with Kia-style event-triggered communication
% ========================================================================
function Xi_final = ip_ac_consensus_kia(Pi, L, Kappa_P, t_step, max_iter, tol, N_degree, epsilon_i)
Xi = Pi;
Xi_hat = Pi;

for iter = 1:max_iter
    Xi_prev = Xi;

    L_Xi_hat = laplacian_multiply_agent_dim_local(Xi_hat, L);
    Xi = Xi - t_step * Kappa_P * L_Xi_hat;

    Xi_hat = kia_trigger_update(Xi, Xi_hat, L, N_degree, epsilon_i, iter);

    if max(abs(Xi(:) - Xi_prev(:))) < tol
        break;
    end
end

Xi_final = Xi;
end

%% ========================================================================
% Kia et al. connected-undirected-graph trigger update
% ========================================================================
function Xi_hat = kia_trigger_update(Xi, Xi_hat, L, N_degree, epsilon_i, iter)
[~, AgentQuantity, NumInducingPoints] = size(Xi);

for agent_i = 1:AgentQuantity
    neighbor_idx = (L(agent_i,:) < 0);
    neighbor_list = find(neighbor_idx);
    d_i = N_degree(agent_i);
    if d_i <= 0
        continue;
    end

    E_i = Xi_hat(:,agent_i,:) - Xi(:,agent_i,:);
    lhs = squeeze(sum(E_i.^2, 1));
    lhs = lhs(:).';

    neighbor_term = zeros(1, NumInducingPoints);
    for jj = 1:numel(neighbor_list)
        nb = neighbor_list(jj);
        D_ij = Xi_hat(:,agent_i,:) - Xi_hat(:,nb,:);
        neighbor_term = neighbor_term + squeeze(sum(D_ij.^2, 1)).';
    end

    rhs = (1/(4*d_i)) * neighbor_term + (1/(4*d_i)) * epsilon_i(agent_i)^2;
    trigger_idx = lhs > rhs;

    % Kia paper assumes t_1^i = 0. In the discrete simulation, this is the
    % initial broadcast of all inducing-point components.
    if iter == 1
        trigger_idx(:) = true;
    end

    if any(trigger_idx)
        Xi_hat(:,agent_i,trigger_idx) = Xi(:,agent_i,trigger_idx);
    end
end
end

%% ========================================================================
% Dataset-style decoder: Xi_final -> phi -> shared MaskedGP
% ========================================================================
function MaskedGP = build_shared_maskedgp_from_xi(Xi_final, method, AgentQuantity, M, p_dim, ...
    InducingPoints_Coordinates, SigmaF, SigmaL, x_dim, y_dim)

method = lower(method);
phi = zeros(y_dim, M);

if p_dim < 2*y_dim
    error('p_dim=%d is too small for y_dim=%d.', p_dim, y_dim);
end

for d = 1:y_dim
    xi1 = squeeze(Xi_final(2*d-1, 1, :))';
    xi2 = squeeze(Xi_final(2*d,   1, :))';

    switch method
        case {'poe','gpoe','bcm','rbcm'}
            den = xi2;
            small_den = abs(den) < 1e-4;
            den(small_den & den >= 0) = 1e-4;
            den(small_den & den <  0) = -1e-4;
            phi(d,:) = xi1 ./ den;

        case 'moe'
            % Kept consistent with the dataset version.
            phi(d,:) = xi1 / AgentQuantity;

        otherwise
            error('Unknown aggregation method: %s.', method);
    end
end

phi(~isfinite(phi)) = 0;

MaskedGP = LocalGP_MultiOutput(x_dim, y_dim, M, 1e-6, SigmaF, SigmaL);
MaskedGP.add_Alldata(InducingPoints_Coordinates, phi);
end

%% ========================================================================
% Apply graph Laplacian along the agent dimension of a 3D tensor
% ========================================================================
function L_X = laplacian_multiply_agent_dim_local(X, L)
[p_dim, agent_quantity, num_points] = size(X);
X_agent_first = permute(X, [2, 1, 3]);
X_flat = reshape(X_agent_first, agent_quantity, []);
L_X_flat = L * X_flat;
L_X_agent_first = reshape(L_X_flat, agent_quantity, p_dim, num_points);
L_X = permute(L_X_agent_first, [2, 1, 3]);
end
