function [TrackingError_vector, t_set] = run_simulation_test_point( ...
    CurrentMode, SaveFolderName, SaveFileName, use_formation, ...
    InitialSeedOverride, OfflineDataQuantityOverride, ...
    UnknownScaleOverride, DisturbanceScaleOverride)
if nargin < 4, use_formation = true; end
if nargin < 5 || isempty(InitialSeedOverride)
    InitialSeedOverride = [];
end
if nargin < 6 || isempty(OfflineDataQuantityOverride)
    OfflineDataQuantityOverride = [];
end
if nargin < 7 || isempty(UnknownScaleOverride)
    UnknownScaleOverride = [];
end
if nargin < 8 || isempty(DisturbanceScaleOverride)
    DisturbanceScaleOverride = [];
end
UnknownScale = 0.2;
DisturbanceScale = 0.1;
if ~isempty(UnknownScaleOverride)
    UnknownScale = UnknownScaleOverride;
end
if ~isempty(DisturbanceScaleOverride)
    DisturbanceScale = DisturbanceScaleOverride;
end
%% 1. System Parameters
SystemOrder = 2; q_dim = 2; x_dim = q_dim * SystemOrder;
m1 = 1; m2 = 1; L1 = 1; L2 = 1; g = 9.8;
AgentQuantity = 6; LeaderQuantity = 1;

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
Lambda = [zeros(SystemOrder-2,1), eye(SystemOrder-2);
          -lambda_set(1)/lambda_n, -lambda_set(2:end-1)/lambda_n];
Qes = eye(SystemOrder-1);
Pes = care(Lambda, [], Qes);
Pe  = kron(Pes, eye(AgentQuantity));
Qe  = kron(Qes, eye(AgentQuantity));
q   = (L + diag(B)) \ ones(AgentQuantity,1);
Pr  = diag(1./q);
Qr  = Pr*(L+diag(B)) + (L+diag(B))'*Pr;
t_vec = [zeros(SystemOrder-2,1); 1/lambda_n];
Phi = Pr*kron(lambda_vector'*t_vec, eye(AgentQuantity));
Psi = Pr*kron(lambda_vector'*Lambda, eye(AgentQuantity)) + ...
      kron(t_vec'*Pes, eye(AgentQuantity));
Pz  = blkdiag(Pr, Pe);
eig_Pz = eig(Pz);
max_eig_Pz = max(eig_Pz);
min_eig_Pz = min(eig_Pz);
Qz  = [c*lambda_n*Qr-2*Phi, -Psi; -Psi', Qe];
eig_Qz = eig(Qz);
min_eig_Qz = min(eig_Qz);
if ~(all(real(eig_Qz)>0) && all(real(eig(Lambda))<0))
    error('Controller is not stable!');
end

xi = 2 * lambda_n / min_eig_Qz * norm(Pr * (L + diag(B)));
chi = sqrt((1 + norm([t_vec, Lambda])^2) * max_eig_Pz / min_eig_Pz) * ...
      norm(inv(L + diag(B)));

%% 4. Time
t_start = 0; t_end = 10; t_step = 0.01;
t_end_override = str2double(getenv('CONTROL_SIM_END_TIME'));
if isfinite(t_end_override) && t_end_override > 0
    t_end = t_end_override;
end
t_set = t_start:t_step:t_end;
T = numel(t_set);
OnlineUpdateInterval = 0.1;
OnlineUpdateStep = max(1,round(OnlineUpdateInterval/t_step));

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
SigmaF = 1; SigmaL = 0.5*ones(x_dim,1);
GP_tau = 1e-8; GP_delta = 0.1; y_dim = q_dim;
DomainScale = 1.5;
MaxDataQuantity_set     = 600*ones(AgentQuantity,1);
DefaultOfflineDataQuantityPerAgent = 350;
SigmaN_set = 0.01*ones(AgentQuantity,1);
prior_var = SigmaF^2;
AggregationParameters = control_aggregation_parameters();

mode_lower = lower(CurrentMode);
dac_methods = {'poe','gpoe','moe','bcm','rbcm'};
ac_methods  = {'poe_ac','gpoe_ac','moe_ac','bcm_ac','rbcm_ac'};
offline_dac_methods = strcat(dac_methods, '_offline');
offline_ac_methods = strcat(ac_methods, '_offline');
is_local_mode = ismember(mode_lower, {'local','local_offline'});
is_offline_mode = ismember(mode_lower, [offline_dac_methods, ...
    offline_ac_methods, {'local_offline'}]);
if is_offline_mode
    if isempty(OfflineDataQuantityOverride)
        OfflineDataQuantity_set = ...
            DefaultOfflineDataQuantityPerAgent*ones(AgentQuantity,1);
    elseif isscalar(OfflineDataQuantityOverride)
        % Match the advisor's manipulator shared code: a scalar is the
        % number of offline samples owned by every agent.
        OfflineDataQuantity_set = max(0,round(OfflineDataQuantityOverride)) ...
            *ones(AgentQuantity,1);
    elseif numel(OfflineDataQuantityOverride) == AgentQuantity
        OfflineDataQuantity_set = OfflineDataQuantityOverride(:);
    else
        error(['OfflineDataQuantityOverride must be a scalar or have ' ...
            'one entry per agent.']);
    end
else
    OfflineDataQuantity_set = zeros(AgentQuantity,1);
end
OfflineDataQuantity = OfflineDataQuantity_set;
OfflineDataQuantityTotal = sum(OfflineDataQuantity_set);

% Keep the offline dataset stream independent of the initial-state stream.
% This makes the dataset vary across Monte Carlo seeds while every method
% within the same seed still starts from exactly the same plant state.
if isempty(InitialSeedOverride)
    OfflineDataSeed = 10042;
else
    OfflineDataSeed = 10000+round(InitialSeedOverride);
end
rng(OfflineDataSeed,'twister');

LocalGP_set = cell(AgentQuantity,1);
for n = 1:AgentQuantity
    LocalGP_set{n} = LocalGP_MultiOutput(x_dim, y_dim, ...
        MaxDataQuantity_set(n), SigmaN_set(n), SigmaF, SigmaL);
    X_in = 2*(rand(x_dim,OfflineDataQuantity_set(n))-0.5)*DomainScale;
    Y_in = UnknownScale * Manipulator_2D_2DoF_UnknownDynamics(X_in);
    Y_in = Y_in + DisturbanceScale * SigmaN_set(n)*randn(size(Y_in));
    LocalGP_set{n}.add_Alldata(X_in, Y_in);
    LocalGP_set{n}.tau = GP_tau;
    LocalGP_set{n}.delta = GP_delta;
end

%% 7. Online Learning ET Setup
Bidirection_NeighbourSet = cell(AgentQuantity, 1);
Sigma_update_aggregation_set = nan(AgentQuantity, 1);

for AgentNr = 1:AgentQuantity
    AgentNeighbourSet = MultiAgentSystem.Agent_Topology.NeighbourSet{AgentNr};

    for NeighbourNr = numel(AgentNeighbourSet):-1:1
        NeighbourAgentNr = AgentNeighbourSet(NeighbourNr);
        if isempty(find(MultiAgentSystem.Agent_Topology.NeighbourSet{NeighbourAgentNr} == AgentNr, 1))
            AgentNeighbourSet(NeighbourNr) = [];
        end
    end

    Bidirection_NeighbourSet{AgentNr} = AgentNeighbourSet;

    Sigma_update_set = nan(numel(AgentNeighbourSet)+1, 1);
    Sigma_update_set(1) = LocalGP_set{AgentNr}.SigmaN;

    for k = 1:numel(AgentNeighbourSet)
        NeighbourAgentNr = AgentNeighbourSet(k);
        Sigma_update_set(k+1) = LocalGP_set{NeighbourAgentNr}.SigmaF;
    end

    Sigma_update_aggregation_set(AgentNr) = ...
        sqrt(1 / (sum(Sigma_update_set.^(-2)) / numel(Sigma_update_set)));
end

beta = 0;
gamma = 0.0005;
for LocalGP_Nr = 1:AgentQuantity
    if LocalGP_set{LocalGP_Nr}.DataQuantity > 0
        [~,~,~,beta_i,~,~] = LocalGP_set{LocalGP_Nr}.predict(zeros(x_dim,1));
        beta = max(beta, beta_i);
    end
end

eta_underline_set = sqrt(beta) * Sigma_update_aggregation_set + gamma;

vartheta_bar = xi * chi * norm( ...
    (eye(AgentQuantity) - diag(B)) * ones(AgentQuantity,1) * Fl ...
    + eta_underline_set);

%% 8. Setup - DAC Zeta 初始化
L_lap   = MultiAgentSystem.Agent_Topology.LaplacianMatrix;
positive_laplacian_eigenvalues = eig(L_lap);
positive_laplacian_eigenvalues = positive_laplacian_eigenvalues( ...
    positive_laplacian_eigenvalues > 1e-10);
TPQueryUpdateInterval = get_env_positive_scalar_tp( ...
    'CONTROL_TP_QUERY_UPDATE_INTERVAL',OnlineUpdateInterval);
TPQueryUpdateStep = max(1,round(TPQueryUpdateInterval/t_step));

% DAC has ten chronological updates over every 0.1 s window. AC instead
% solves each new static reference with ten fixed inner iterations.
TPDACFixedRounds = get_env_positive_integer_tp( ...
    'CONTROL_TP_DAC_FIXED_ROUNDS',TPQueryUpdateStep);
TPACFixedRounds = get_env_positive_integer_tp( ...
    'CONTROL_TP_AC_FIXED_ROUNDS',TPQueryUpdateStep);
if mod(TPDACFixedRounds,TPQueryUpdateStep) ~= 0
    error(['TP-DAC fixed rounds must be an integer multiple of the number of ' ...
        'physical steps in one query window (%d).'],TPQueryUpdateStep);
end
TPDACStepsPerPhysicalStep = TPDACFixedRounds/TPQueryUpdateStep;
TPACRoundsPerReferenceUpdate = TPACFixedRounds;
TPACRoundsPerPhysicalStep = NaN; % legacy field; AC uses an inner solve

% DAC retains its paper-style gain. AC has a separate discrete step so its
% ten-round comparison can use the same value as IP-AC.
default_tp_consensus_step = 1.9/(min(positive_laplacian_eigenvalues) ...
    + max(positive_laplacian_eigenvalues));
TPConsensusKappaP = get_env_positive_scalar_tp( ...
    'CONTROL_TP_CONSENSUS_KAPPA_P',default_tp_consensus_step/t_step);
TPDACConsensusStep = t_step*TPConsensusKappaP;
TPACConsensusStep = get_env_positive_scalar_tp( ...
    'CONTROL_TP_AC_CONSENSUS_STEP',0.2);
if TPDACConsensusStep >= 2/max(positive_laplacian_eigenvalues)
    error(['Unstable TP consensus gain: t_step*Kappa_P must be smaller ' ...
        'than 2/lambda_max(L).']);
end
if TPACConsensusStep >= 2/max(positive_laplacian_eigenvalues)
    error(['Unstable TP-AC consensus gain: its discrete step must be ' ...
        'smaller than 2/lambda_max(L).']);
end

% Compatibility fields retained in saved MAT files. Neither TP method
% uses a tolerance or an extra initialization burst.
TPDACInitialRounds = 0;
ACMaxIterations = TPACFixedRounds;
ACConsensusTolerance = NaN;
TPACIterationPolicy = 'fixed-reference-inner-rounds';
TPConsensusTiming = 'dac-continuous-ac-inner-fixed';
TPImplementationVersion = 3;

p_dim_tp = 2 * y_dim;
% TP uses one information vector per agent:
%   column i = information produced by GP_i at agent i's own state x_i.
% No agent evaluates its GP at a neighbour's query point.
% Hence the consensus state has size p_dim_tp x AgentQuantity only.
TPQueryPolicy = 'own-agent-query';
TPConsensusStateSize = [p_dim_tp,AgentQuantity];
Zeta_vector = zeros(p_dim_tp,AgentQuantity);
P_tp_snapshot = zeros(p_dim_tp,AgentQuantity);
Xi_tp_ac_state = zeros(p_dim_tp,AgentQuantity);
tp_snapshot_initialized = false;
tp_ac_round_in_snapshot = 0;

%% 8. Initial State
if ~isempty(InitialSeedOverride)
    rng(InitialSeedOverride);
else
    rng(42);  % 固定初始状态seed，确保所有方法一致
end
x_all = rand(x_dim*AgentQuantity, 1);
x_all_set = nan(x_dim*AgentQuantity, T);
x_all_set(:,1) = x_all;
%fprintf('[%s] 初始状态 x_all(1:4): %.4f %.4f %.4f %.4f\n', CurrentMode, x_all(1:4));
vartheta_all_set = nan(x_dim*AgentQuantity, T);
vartheta_all_set(:,1) = x_all - s_all_set(:,1) - kron(ones(AgentQuantity,1), xl_set(:,1));

f_hat_matrix  = zeros(y_dim, AgentQuantity);
f_true_matrix = zeros(y_dim, AgentQuantity);
TrackingError_vector = zeros(1, T);
f_hat_all_set  = nan(y_dim, AgentQuantity, T);
f_true_all_set = nan(y_dim, AgentQuantity, T);
prediction_error_norm_vector = nan(1,T);
tp_consensus_error_vector = nan(1,T);
tp_consensus_round_count_set = zeros(1,T);
tp_broadcast_count_per_agent = zeros(AgentQuantity,1);
tp_ac_round_error_set = nan(TPACFixedRounds+1,T);

online_trigger_set = zeros(AgentQuantity, T);
online_trigger_count = zeros(AgentQuantity, 1);

%% 9. Control Loop
opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-3, ...
    'InitialStep', t_step, 'MaxStep', t_step, 'Refine', 1);
tic;
for t_Nr = 1:T-1
    t = t_set(t_Nr);
    x_l_r        = xlr_set(:, t_Nr);
    x_all        = x_all_set(:, t_Nr);
    x_all_matrix = reshape(x_all, x_dim, AgentQuantity);
    x_all_cell   = ET_MAS_GP_Leader_vector2cell(x_all, AgentQuantity, 1);
    s_all        = s_all_set(:, t_Nr);
    s_r_all      = sr_all_set(:, t_Nr);
    s_r_cell     = ET_MAS_GP_Leader_vector2cell(s_r_all, AgentQuantity, 1);
    x_tilde_all  = x_all - s_all;
    x_tilde_cell = ET_MAS_GP_Leader_vector2cell(x_tilde_all, AgentQuantity, SystemOrder);
    vartheta_all  = vartheta_all_set(:, t_Nr);
    vartheta_cell = ET_MAS_GP_Leader_vector2cell(vartheta_all, AgentQuantity, SystemOrder);

    TrackingError_vector(t_Nr) = norm(vartheta_all);

    % Match the inducing-point baseline: the current measurement is added
    % before building P and before computing the current control input.
    is_online_update = mod(t_Nr-1,OnlineUpdateStep) == 0;
    if is_online_update && ~strcmpi(CurrentMode, 'exact') && ~is_offline_mode
        [LocalGP_set, online_trigger_set, online_trigger_count] = ...
            apply_time_triggered_online_learning(LocalGP_set, ...
            online_trigger_set, online_trigger_count, t_Nr, x_all_matrix, ...
            UnknownScale, DisturbanceScale);
    end

    [phi_cell, r_matrix, e_cell] = Manipulator_2D_2DoF_ConsensusLaw( ...
        vartheta_cell, x_tilde_cell, x_l_r, MultiAgentSystem, c, lambda_set, s_r_cell);

    AgentState_matrix = x_all_matrix;
    is_tp_dac_mode = ismember(mode_lower,[dac_methods,offline_dac_methods]);
    is_tp_ac_mode = ismember(mode_lower,[ac_methods,offline_ac_methods]);
    refresh_tp_snapshot = is_tp_ac_mode && ...
        (~tp_snapshot_initialized || mod(t_Nr-1,TPQueryUpdateStep) == 0);
    if refresh_tp_snapshot
        tp_base = strrep(strrep(mode_lower,'_offline',''),'_ac','');
        P_tp_snapshot = build_tp_own_query_prediction_info( ...
            LocalGP_set,AgentState_matrix,AgentQuantity,y_dim, ...
            tp_base,prior_var,AggregationParameters);
        tp_snapshot_initialized = true;
        Xi_tp_ac_state = P_tp_snapshot;
        tp_ac_round_in_snapshot = 0;
    end

    if is_tp_dac_mode
        % TP-DAC follows (25)/(29) of the distributed-GP paper:
        %   1) every physical step, GP_i is queried only at x_i[k];
        %   2) Zeta is retained over the complete trajectory;
        %   3) one chronological communication update is performed here.
        % There is no snapshot reset, tolerance, or terminal condition.
        base = strrep(mode_lower, '_offline', '');
        P_tp_current = build_tp_own_query_prediction_info( ...
            LocalGP_set,AgentState_matrix,AgentQuantity,y_dim, ...
            base,prior_var,AggregationParameters);
        rounds_this_step = TPDACStepsPerPhysicalStep;
        for dac_round = 1:rounds_this_step
            Xi_tp = P_tp_current-Zeta_vector;
            Zeta_vector = Zeta_vector + TPDACConsensusStep* ...
                laplacian_multiply_tp(Xi_tp,L_lap);
        end

        Xi_tp = P_tp_current-Zeta_vector;
        P_average = mean(P_tp_current,2);
        tp_consensus_error_vector(t_Nr) = ...
            max(abs(Xi_tp-P_average),[],'all');
        for agent_i = 1:AgentQuantity
            f_hat_matrix(:,agent_i) = decode_tp_prediction_info( ...
                Xi_tp(:,agent_i),base,AgentQuantity,y_dim);
            f_hat_matrix(:,agent_i) = max(-30,min(30, ...
                f_hat_matrix(:,agent_i)));
        end
        % One round means one p-vector broadcast by each agent.
        tp_consensus_round_count_set(t_Nr) = rounds_this_step;
        tp_broadcast_count_per_agent = tp_broadcast_count_per_agent ...
            + rounds_this_step;

    elseif is_tp_ac_mode
        % TP-AC is the static-average comparison baseline. P is sampled
        % and Xi is reset once per 0.1 s window. The new reference is held
        % fixed while exactly ten inner AC iterations are completed.
        base = strrep(strrep(mode_lower,'_offline',''),'_ac','');
        if refresh_tp_snapshot
            step_round_count = TPACFixedRounds;
        else
            step_round_count = 0;
        end

        P_average = mean(P_tp_snapshot,2);
        if refresh_tp_snapshot
            initial_error = Xi_tp_ac_state-P_average;
            tp_ac_round_error_set(1,t_Nr) = ...
                sqrt(mean(initial_error.^2,'all'));
        end

        for ac_round = 1:step_round_count
            Xi_tp_ac_state = Xi_tp_ac_state-TPACConsensusStep* ...
                laplacian_multiply_tp(Xi_tp_ac_state,L_lap);
            tp_ac_round_in_snapshot = tp_ac_round_in_snapshot+1;
            round_error = Xi_tp_ac_state-P_average;
            tp_ac_round_error_set(ac_round+1,t_Nr) = ...
                sqrt(mean(round_error.^2,'all'));
        end

        tp_consensus_error_vector(t_Nr) = max(abs( ...
            Xi_tp_ac_state-P_average),[],'all');
        for agent_i = 1:AgentQuantity
            f_hat_matrix(:,agent_i) = decode_tp_prediction_info( ...
                Xi_tp_ac_state(:,agent_i),base,AgentQuantity,y_dim);
            f_hat_matrix(:,agent_i) = max(-30,min(30, ...
                f_hat_matrix(:,agent_i)));
        end
        tp_consensus_round_count_set(t_Nr) = step_round_count;
        tp_broadcast_count_per_agent = tp_broadcast_count_per_agent ...
            + step_round_count;

    elseif is_local_mode
        for n = 1:AgentQuantity
            [mu_n,~] = predict_gp_mean_variance( ...
                LocalGP_set{n},AgentState_matrix(:,n));
            mu_n = max(-30, min(30, mu_n));
            f_hat_matrix(:,n) = mu_n;
        end

    elseif strcmpi(CurrentMode,'exact')
        for n = 1:AgentQuantity
            f_hat_matrix(:,n) = UnknownScale * ...
                Manipulator_2D_2DoF_UnknownDynamics(AgentState_matrix(:,n));
        end
    end

    for n = 1:AgentQuantity
        f_true_matrix(:,n) = UnknownScale * ...
            Manipulator_2D_2DoF_UnknownDynamics(x_all_matrix(:,n));
    end
    f_hat_all_set(:,:,t_Nr)  = f_hat_matrix;
    f_true_all_set(:,:,t_Nr) = f_true_matrix;
    prediction_error_norm_vector(t_Nr) = ...
        norm(f_true_matrix(:)-f_hat_matrix(:));

    u_cell = Manipulator_2D_2DoF_get_u_cell(x_all_cell, phi_cell, f_hat_matrix, L1, L2, m1, m2);

    [~, x_all_temp] = ode45( ...
        @(t,x) Manipulator_2D_2DoF_MultiAgent_DynamicFunction( ...
        t, x, u_cell, L1, L2, m1, m2, UnknownScale), ...
        [t, t+t_step], x_all, opts);
    x_all_next = x_all_temp(end,:)';
    x_all_set(:, t_Nr+1) = x_all_next;
    vartheta_all_set(:, t_Nr+1) = x_all_next - s_all_set(:,t_Nr+1) - ...
        kron(ones(AgentQuantity,1), xl_set(:,t_Nr+1));

    % fprintf('t = %6.4f\n', t);
end
TrackingError_vector(end) = norm(vartheta_all_set(:,end));
evaluation_mask_after10s = t_set >= 10;
if any(evaluation_mask_after10s)
    MaxTrackingErrorAfter10s = max( ...
        TrackingError_vector(evaluation_mask_after10s),[],'omitnan');
else
    MaxTrackingErrorAfter10s = NaN;
end

x_all_matrix_end = reshape(x_all_set(:,end), x_dim, AgentQuantity);
for n = 1:AgentQuantity
    f_true_all_set(:,n,end) = UnknownScale * ...
        Manipulator_2D_2DoF_UnknownDynamics(x_all_matrix_end(:,n));
    f_hat_all_set(:,n,end)  = f_hat_matrix(:,n);
end
prediction_error_norm_vector(end) = norm(reshape( ...
    f_true_all_set(:,:,end)-f_hat_all_set(:,:,end),[],1));

elapsed_time = toc;

fprintf('==================================================');
fprintf('Mode: %s', CurrentMode);
fprintf('Formation: %d', use_formation);
fprintf('Total simulation time: %.2f s', elapsed_time);
fprintf('Max tracking error for t >= 10 s: %.6f', ...
    MaxTrackingErrorAfter10s);
fprintf('Online learning time trigger:');
fprintf('  Total triggers: %d', sum(online_trigger_count));
fprintf('  Average triggers: %.2f / agent', mean(online_trigger_count));
fprintf('==================================================');

%% 10. Save
if nargin >= 3 && ~isempty(SaveFolderName) && ~isempty(SaveFileName)
    if ~exist(SaveFolderName,'dir'), mkdir(SaveFolderName); end
    save(fullfile(SaveFolderName,[SaveFileName,'.mat']), ...
        't_set','TrackingError_vector','MaxTrackingErrorAfter10s', ...
        'CurrentMode','use_formation',...
        'f_hat_all_set','f_true_all_set','prediction_error_norm_vector', ...
        'vartheta_all_set', ...
        'online_trigger_set','online_trigger_count','eta_underline_set', ...
        'vartheta_bar','elapsed_time','UnknownScale','DisturbanceScale', ...
        'OfflineDataQuantity','OfflineDataQuantity_set', ...
        'OfflineDataQuantityTotal','is_offline_mode', ...
        'OfflineDataSeed', ...
        'OnlineUpdateInterval','OnlineUpdateStep', ...
        'AggregationParameters', ...
        'TPDACStepsPerPhysicalStep','TPDACFixedRounds', ...
        'TPDACConsensusStep','TPConsensusKappaP', ...
        'TPDACInitialRounds','TPACConsensusStep','ACMaxIterations', ...
        'TPACIterationPolicy','TPACFixedRounds', ...
        'TPQueryUpdateInterval','TPQueryUpdateStep', ...
        'TPQueryPolicy','TPConsensusStateSize', ...
        'TPACRoundsPerPhysicalStep','TPACRoundsPerReferenceUpdate', ...
        'TPConsensusTiming', ...
        'TPImplementationVersion', ...
        'ACConsensusTolerance','tp_consensus_error_vector', ...
        'tp_consensus_round_count_set','tp_broadcast_count_per_agent', ...
        'tp_ac_round_error_set');
end
end

function L_Xi = laplacian_multiply_tp(Xi,L)
% Xi columns are model/communication agents. Column i of the result is
% sum_r L(i,r)*Xi_r, with the same column layout as Xi.
L_Xi = Xi*L.';
end

function value = get_env_positive_integer_tp(name,default_value)
value = get_env_positive_scalar_tp(name,default_value);
value = round(value);
end

function value = get_env_positive_scalar_tp(name,default_value)
value = str2double(getenv(name));
if ~(isfinite(value) && value > 0)
    value = default_value;
end
end

function P_tp = build_tp_own_query_prediction_info( ...
    LocalGP_set,AgentState_matrix,AgentQuantity,y_dim,method,prior_var, ...
    aggregation_cfg)
% Each agent evaluates only its own GP at its own state:
%       P_i[k] = information(GP_i,x_i[k]).
% The resulting p-by-N matrix is the single TP consensus problem.

P_tp = zeros(2*y_dim, AgentQuantity);

for agent_i = 1:AgentQuantity
    [mu_i,var_i] = predict_gp_mean_variance( ...
        LocalGP_set{agent_i},AgentState_matrix(:,agent_i));

    mu_i = mu_i(:);
    if isscalar(var_i)
        var_i = var_i*ones(y_dim,1);
    else
        var_i = var_i(:);
    end

    for d = 1:y_dim
        sv = max(var_i(d),aggregation_cfg.posterior_var_floor);
        raw_beta = 0.5*(log(prior_var)-log(sv));
        beta_gpoe = max(min(raw_beta,aggregation_cfg.gpoe_beta_max),eps);
        beta_rbcm = max(min(raw_beta,aggregation_cfg.rbcm_beta_max),eps);

        switch lower(method)
            case 'moe'
                P_tp(2*d-1,agent_i) = AgentQuantity*mu_i(d);
                P_tp(2*d,  agent_i) = AgentQuantity*(sv+mu_i(d)^2);

            case 'gpoe'
                P_tp(2*d-1,agent_i) = AgentQuantity*beta_gpoe*mu_i(d)/sv;
                P_tp(2*d,  agent_i) = AgentQuantity*beta_gpoe/sv;

            case 'poe'
                P_tp(2*d-1,agent_i) = AgentQuantity*mu_i(d)/sv;
                P_tp(2*d,  agent_i) = AgentQuantity/sv;

            case 'bcm'
                P_tp(2*d-1,agent_i) = AgentQuantity*mu_i(d)/sv;
                P_tp(2*d,  agent_i) = AgentQuantity/sv- ...
                    aggregation_cfg.bcm_prior_scale* ...
                    (AgentQuantity-1)/prior_var;

            case 'rbcm'
                P_tp(2*d-1,agent_i) = ...
                    AgentQuantity*beta_rbcm*mu_i(d)/sv;
                P_tp(2*d,agent_i) = AgentQuantity*beta_rbcm/sv+ ...
                    (1-AgentQuantity*beta_rbcm)/prior_var;

            otherwise
                error('Unknown TP aggregation method: %s', method);
        end
    end
end
end

function mu_hat = decode_tp_prediction_info(p_vec, method, AgentQuantity, y_dim)

mu_hat = zeros(y_dim,1);

for d = 1:y_dim
    xi1 = p_vec(2*d-1);
    xi2 = p_vec(2*d);

    if ismember(lower(method), {'gpoe','poe','bcm','rbcm'})
        mu_hat(d) = xi1 / max(abs(xi2), 1e-4);
    else
        mu_hat(d) = xi1 / AgentQuantity;
    end
end

mu_hat(~isfinite(mu_hat)) = 0;

end


function [LocalGP_set, online_trigger_set, online_trigger_count] = ...
    apply_online_learning_et( ...
    LocalGP_set, online_trigger_set, online_trigger_count, ...
    t_Nr, x_all_matrix, r_matrix, e_cell, ...
    AgentQuantity, y_dim, beta, gamma, ...
    Bidirection_NeighbourSet, eta_underline_set, ...
    MultiAgentSystem, Fl, xi, chi, vartheta_bar)

mu_cell = cell(AgentQuantity, AgentQuantity);
var_matrix = nan(AgentQuantity, AgentQuantity);
eta_aggregated_vector = nan(AgentQuantity, 1);

% Compute current aggregated eta for each controlled agent.
% Formation convention:
%   agent i needs f(x_i), therefore all GP models are queried at x_i.
for AgentNr = 1:AgentQuantity

    x_i = x_all_matrix(:, AgentNr);

    [mu_cell{AgentNr,AgentNr}, var_matrix(AgentNr,AgentNr), ~] = ...
        Manipulator_2D_2DoF_LocalPrediction( ...
        x_i, AgentNr, LocalGP_set, beta, gamma, y_dim);

    AgentBidirection_NeighbourSet = Bidirection_NeighbourSet{AgentNr};

    for k = 1:numel(AgentBidirection_NeighbourSet)
        NeighbourAgentNr = AgentBidirection_NeighbourSet(k);

        % Use neighbour's GP model, but evaluate at agent i's state x_i.
        [mu_cell{AgentNr,NeighbourAgentNr}, ...
         var_matrix(AgentNr,NeighbourAgentNr), ~] = ...
            Manipulator_2D_2DoF_LocalPrediction( ...
            x_i, NeighbourAgentNr, LocalGP_set, beta, gamma, y_dim);
    end

    [~, eta_aggregated_i] = ...
        ET_MAS_GP_Leader_GPAggregation_SingleAgent( ...
        AgentNr, AgentBidirection_NeighbourSet, ...
        var_matrix(AgentNr,:), mu_cell(AgentNr,:), beta, gamma);

    eta_aggregated_vector(AgentNr) = eta_aggregated_i;
end

% Online data-trigger decision.
for AgentNr = 1:AgentQuantity
    online_trigger_set(AgentNr, t_Nr) = ...
        Manipulator_2D_2DoF_DistributedET( ...
        AgentNr, r_matrix, e_cell, ...
        eta_underline_set, eta_aggregated_vector, ...
        MultiAgentSystem, Fl, xi, chi, vartheta_bar);
end

% Add online data to triggered local GPs.
for AgentNr = 1:AgentQuantity
    if online_trigger_set(AgentNr, t_Nr) == 1

        x_i = x_all_matrix(:, AgentNr);

        y_i = Manipulator_2D_2DoF_UnknownDynamics(x_i) + ...
              LocalGP_set{AgentNr}.SigmaN * randn(y_dim, 1);

        if LocalGP_set{AgentNr}.DataQuantity >= LocalGP_set{AgentNr}.MaxDataQuantity
            LocalGP_set{AgentNr}.downdateParam(1);
        end

        LocalGP_set{AgentNr}.addPoint(x_i, y_i);

        online_trigger_count(AgentNr) = online_trigger_count(AgentNr) + 1;
    end
end

end
