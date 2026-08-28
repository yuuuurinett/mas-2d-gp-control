function [TrackingError_vector, t_set] = run_simulation_nbr( ...
    CurrentMode,SaveFolderName,SaveFileName,use_formation, ...
    InitialSeedOverride,OfflineDataQuantityOverride, ...
    UnknownScaleOverride,DisturbanceScaleOverride)
if nargin < 4, use_formation = true; end
if nargin < 5 || isempty(InitialSeedOverride)
    InitialSeedOverride = [];
end
if nargin < 6 || isempty(OfflineDataQuantityOverride)
    OfflineDataQuantityOverride = [];
end
if nargin < 7, UnknownScaleOverride = []; end
if nargin < 8, DisturbanceScaleOverride = []; end
rng(0);
UnknownScale = 0.2;
DisturbanceScale = 0.1;
if ~isempty(UnknownScaleOverride), UnknownScale = UnknownScaleOverride; end
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
DefaultOfflineDataQuantity = 350;
SigmaN_set = 0.01*ones(AgentQuantity,1);
prior_var = SigmaF^2;
aggregation_cfg = control_aggregation_parameters();

mode_lower = lower(CurrentMode);
dac_methods = {'poe','gpoe','moe','bcm','rbcm'};
offline_methods = strcat(dac_methods, '_offline');
is_offline_mode = ismember(mode_lower, offline_methods);
base_mode = strrep(mode_lower, '_offline', '');
if is_offline_mode
    if isempty(OfflineDataQuantityOverride)
        OfflineDataQuantity_set = ...
            DefaultOfflineDataQuantity*ones(AgentQuantity,1);
    elseif isscalar(OfflineDataQuantityOverride)
        OfflineDataQuantity_set = ...
            OfflineDataQuantityOverride*ones(AgentQuantity,1);
    else
        OfflineDataQuantity_set = OfflineDataQuantityOverride(:);
        if numel(OfflineDataQuantity_set) ~= AgentQuantity
            error('Offline sample allocation must contain one value per agent.');
        end
    end
else
    OfflineDataQuantity_set = zeros(AgentQuantity,1);
end
OfflineDataQuantity = mean(OfflineDataQuantity_set);

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

% NBR uses self + bidirectional neighbours.
NeighborSet = cell(AgentQuantity,1);
for n = 1:AgentQuantity
    NeighborSet{n} = [n, Bidirection_NeighbourSet{n}(:)'];
end

%% 8. Initial State
if ~isempty(InitialSeedOverride)
    rng(InitialSeedOverride);
else
    rng(42);  % 固定初始状态seed，确保所有方法一致
end
x_all = rand(x_dim*AgentQuantity, 1);
x_all_set = nan(x_dim*AgentQuantity, T);
x_all_set(:,1) = x_all;
vartheta_all_set = nan(x_dim*AgentQuantity, T);
vartheta_all_set(:,1) = x_all - s_all_set(:,1) - kron(ones(AgentQuantity,1), xl_set(:,1));

f_hat_matrix  = zeros(y_dim, AgentQuantity);
f_true_matrix = zeros(y_dim, AgentQuantity);
TrackingError_vector = zeros(1, T);
f_hat_all_set  = nan(y_dim, AgentQuantity, T);
f_true_all_set = nan(y_dim, AgentQuantity, T);
prediction_error_norm_vector = nan(1,T);

online_trigger_set = zeros(AgentQuantity, T);
online_trigger_count = zeros(AgentQuantity, 1);

%% 8. Control Loop
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

    % Add the current measurement before the current GP prediction.
    is_online_update = mod(t_Nr-1,OnlineUpdateStep) == 0;
    if is_online_update && ~is_offline_mode && ~strcmpi(CurrentMode,'exact')
        [LocalGP_set, online_trigger_set, online_trigger_count] = ...
            apply_time_triggered_online_learning(LocalGP_set, ...
            online_trigger_set, online_trigger_count, t_Nr, x_all_matrix, ...
            UnknownScale, DisturbanceScale);
    end

    [phi_cell, r_matrix, e_cell] = Manipulator_2D_2DoF_ConsensusLaw( ...
        vartheta_cell, x_tilde_cell, x_l_r, MultiAgentSystem, c, lambda_set, s_r_cell);

    for n = 1:AgentQuantity
        x_n   = x_all_matrix(:, n);
        nbrs  = NeighborSet{n};
        N_nbr = numel(nbrs);
        mu_nbr  = zeros(y_dim, N_nbr);
        var_nbr = zeros(y_dim, N_nbr);
        for ki = 1:N_nbr
            k = nbrs(ki);
            [mu_k, var_k] = predict_gp_mean_variance(LocalGP_set{k},x_n);
            mu_nbr(:,ki)  = mu_k;
            var_nbr(:,ki) = var_k;
        end
        mu_agg = zeros(y_dim, 1);
        for d = 1:y_dim
            mu_a  = mu_nbr(d,:)';
            var_a = var_nbr(d,:)';
            var_a = max(var_a,aggregation_cfg.posterior_var_floor);
            raw_beta = 0.5*(log(prior_var)-log(var_a));
            beta_gpoe = max(min(raw_beta, ...
                aggregation_cfg.gpoe_beta_max),eps);
            beta_rbcm = max(min(raw_beta, ...
                aggregation_cfg.rbcm_beta_max),eps);
            switch base_mode
                case 'moe'
                    mu_agg(d) = mean(mu_a);
                case 'gpoe'
                    prec = sum(beta_gpoe./var_a);
                    mu_agg(d) = sum(beta_gpoe.*mu_a./var_a) / ...
                        max(prec,aggregation_cfg.precision_floor);
                case 'poe'
                    prec = sum(1./var_a);
                    mu_agg(d) = sum(mu_a./var_a) / max(prec, eps);
                case 'bcm'
                    prec = sum(1./var_a)-aggregation_cfg.bcm_prior_scale* ...
                        (N_nbr-1)/prior_var;
                    mu_agg(d) = sum(mu_a./var_a) / ...
                        max(prec,aggregation_cfg.precision_floor);
                case 'rbcm'
                    prec = sum(beta_rbcm./var_a)+ ...
                        (1-sum(beta_rbcm))/prior_var;
                    mu_agg(d) = sum(beta_rbcm.*mu_a./var_a) / ...
                        max(prec,aggregation_cfg.precision_floor);
            end
        end
        mu_agg = max(-30, min(30, mu_agg));
        f_hat_matrix(:, n) = mu_agg;
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
fprintf('Final tracking error: %.6f', TrackingError_vector(end));
fprintf('Online learning time trigger:');
fprintf('  Total triggers: %d', sum(online_trigger_count));
fprintf('  Average triggers: %.2f / agent', mean(online_trigger_count));
fprintf('==================================================');

%% 9. Save
if nargin >= 3 && ~isempty(SaveFolderName) && ~isempty(SaveFileName)
    if ~exist(SaveFolderName,'dir'), mkdir(SaveFolderName); end
    save(fullfile(SaveFolderName,[SaveFileName,'.mat']), ...
        't_set','TrackingError_vector','CurrentMode','use_formation',...
        'f_hat_all_set','f_true_all_set','prediction_error_norm_vector', ...
        'vartheta_all_set', ...
        'online_trigger_set','online_trigger_count','eta_underline_set', ...
        'vartheta_bar','elapsed_time','UnknownScale','DisturbanceScale', ...
        'aggregation_cfg', ...
        'OfflineDataQuantity','OfflineDataQuantity_set','is_offline_mode', ...
        'OnlineUpdateInterval','OnlineUpdateStep');
end
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
