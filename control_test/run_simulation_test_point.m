function [TrackingError_vector, t_set] = run_simulation_test_point(CurrentMode, SaveFolderName, SaveFileName, use_formation)
if nargin < 4, use_formation = true; end
rng(0);
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
SigmaF = 1; SigmaL = 0.5*ones(x_dim,1);
GP_tau = 1e-8; GP_delta = 0.1; y_dim = q_dim;
DomainScale = 1.5;
MaxDataQuantity_set     = 600*ones(AgentQuantity,1);
OfflineDataQuantity_set = 200*ones(AgentQuantity,1);
SigmaN_set = 0.01*ones(AgentQuantity,1);
prior_var = SigmaF^2;

LocalGP_set = cell(AgentQuantity,1);
for n = 1:AgentQuantity
    LocalGP_set{n} = LocalGP_MultiOutput(x_dim, y_dim, ...
        MaxDataQuantity_set(n), SigmaN_set(n), SigmaF, SigmaL);
    X_in = 2*(rand(x_dim,OfflineDataQuantity_set(n))-0.5)*DomainScale;
    Y_in = Manipulator_2D_2DoF_UnknownDynamics(X_in);
    Y_in = Y_in + SigmaN_set(n)*randn(size(Y_in));
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
    [~,~,~,beta_i,~,~] = LocalGP_set{LocalGP_Nr}.predict(zeros(x_dim,1));
    beta = max(beta, beta_i);
end

eta_underline_set = sqrt(beta) * Sigma_update_aggregation_set + gamma;

vartheta_bar = xi * chi * norm( ...
    (eye(AgentQuantity) - diag(B)) * ones(AgentQuantity,1) * Fl ...
    + eta_underline_set);

%% 8. Setup - DAC Zeta 初始化
dac_methods = {'poe','gpoe','moe','bcm','rbcm'};
ac_methods  = {'poe_ac','gpoe_ac','moe_ac','bcm_ac','rbcm_ac'};
Kappa_P = 100;
L_lap   = MultiAgentSystem.Agent_Topology.LaplacianMatrix;

p_dim_tp = 2 * y_dim;
Zeta_vector = zeros(p_dim_tp, AgentQuantity, AgentQuantity);  % p_dim x model-agent x query-agent

%% 8. Initial State
rng(42);  % 固定初始状态seed，确保所有方法一致
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

online_trigger_set = zeros(AgentQuantity, T);
online_trigger_count = zeros(AgentQuantity, 1);

%% 9. Control Loop
opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
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

    [phi_cell, r_matrix, e_cell] = Manipulator_2D_2DoF_ConsensusLaw( ...
        vartheta_cell, x_tilde_cell, x_l_r, MultiAgentSystem, c, lambda_set, s_r_cell);

    AgentState_matrix = x_all_matrix;

    if ismember(lower(CurrentMode), dac_methods)
        % TP-DAC:
        % For each controlled agent n, build prediction information at x_n.
        % Every GP model k is queried at x_n, not at its own x_k.
        base = lower(CurrentMode);

        for n = 1:AgentQuantity
            x_n = AgentState_matrix(:, n);

            P_tp = build_tp_prediction_info( ...
                x_n, LocalGP_set, AgentQuantity, y_dim, ...
                base, prior_var, SigmaF);

            Xi_tp = P_tp - Zeta_vector(:,:,n);

            L_Xi = zeros(size(Xi_tp));
            for ai = 1:AgentQuantity
                L_Xi(:,ai) = sum(Xi_tp .* reshape(L_lap(ai,:),1,AgentQuantity), 2);
            end

            Zeta_vector(:,:,n) = Zeta_vector(:,:,n) + ...
                t_step * Kappa_P * L_Xi;

            Xi_tp = P_tp - Zeta_vector(:,:,n);

            f_hat_matrix(:,n) = decode_tp_prediction_info( ...
                Xi_tp(:,n), base, AgentQuantity, y_dim);

            f_hat_matrix(:,n) = max(-30, min(30, f_hat_matrix(:,n)));
        end

    elseif ismember(lower(CurrentMode), ac_methods)
        % TP-AC：在当前状态点上构建P矩阵，做AC迭代共识
        % 和dataset的TP-AC一致：Xi初始值=Pi，纯共识迭代收敛到Pi的平均
        base = strrep(lower(CurrentMode), '_ac', '');
        p_dim_tp = 2 * y_dim;
        dac_step_size = 0.01; dac_gain = 10; max_iters_tp = 3000;

        for n = 1:AgentQuantity
            x_n = AgentState_matrix(:, n);

            % 构建P矩阵（当前状态点）
            P_tp = zeros(p_dim_tp, AgentQuantity);
            for k = 1:AgentQuantity
                [mu_k, var_k] = LocalGP_set{k}.predict(x_n);
                for d = 1:y_dim
                    sv   = max(var_k(d), 1e-6);
                    beta = max(min(0.5*(log(prior_var)-log(sv)), 10), eps);
                    switch base
                        case 'moe'
                            P_tp(2*d-1,k) = AgentQuantity * mu_k(d);
                            P_tp(2*d,  k) = AgentQuantity * (sv + mu_k(d)^2);
                        case 'gpoe'
                            P_tp(2*d-1,k) = AgentQuantity * beta * mu_k(d) / sv;
                            P_tp(2*d,  k) = AgentQuantity * beta / sv;
                        case 'poe'
                            P_tp(2*d-1,k) = AgentQuantity * mu_k(d) / sv;
                            P_tp(2*d,  k) = AgentQuantity / sv;
                        case 'bcm'
                            P_tp(2*d-1,k) = AgentQuantity * mu_k(d) / sv;
                            P_tp(2*d,  k) = AgentQuantity / sv - (AgentQuantity-1)/prior_var;
                        case 'rbcm'
                            P_tp(2*d-1,k) = AgentQuantity * beta * mu_k(d) / sv;
                            P_tp(2*d,  k) = AgentQuantity * beta / sv + (1-AgentQuantity*beta)/prior_var;
                    end
                end
            end

            % AC迭代共识：Xi初始值=Pi，收敛到Pi的平均
            Xi_tp = P_tp; dac_iter_tp = 0;
            while dac_iter_tp < max_iters_tp
                dac_iter_tp = dac_iter_tp + 1;
                Xi_tp_prev = Xi_tp;
                L_Xi = zeros(size(Xi_tp));
                for ai = 1:AgentQuantity
                    L_Xi(:,ai) = sum(Xi_tp .* reshape(L_lap(ai,:),1,AgentQuantity), 2);
                end
                Xi_tp = Xi_tp - dac_step_size * dac_gain * L_Xi;
                Xi_mean = mean(Xi_tp, 2);
                if max(abs(Xi_tp - Xi_mean), [], 'all') < 1e-5, break; end
            end

            % 提取agent n的预测
            for d = 1:y_dim
                xi1 = Xi_tp(2*d-1, n);
                xi2 = Xi_tp(2*d,   n);
                if ismember(base, {'gpoe','poe','bcm','rbcm'})
                    f_hat_matrix(d,n) = xi1 / max(abs(xi2), 1e-4);
                else
                    f_hat_matrix(d,n) = xi1 / AgentQuantity;
                end
            end
            f_hat_matrix(:,n) = max(-30, min(30, f_hat_matrix(:,n)));
        end

    elseif strcmpi(CurrentMode,'local')
        for n = 1:AgentQuantity
            [mu_n,~] = LocalGP_set{n}.predict(AgentState_matrix(:,n));
            mu_n = max(-30, min(30, mu_n));
            f_hat_matrix(:,n) = mu_n;
        end

    elseif strcmpi(CurrentMode,'exact')
        for n = 1:AgentQuantity
            f_hat_matrix(:,n) = Manipulator_2D_2DoF_UnknownDynamics(AgentState_matrix(:,n));
        end
    end

    for n = 1:AgentQuantity
        f_true_matrix(:,n) = Manipulator_2D_2DoF_UnknownDynamics(x_all_matrix(:,n));
    end
    f_hat_all_set(:,:,t_Nr)  = f_hat_matrix;
    f_true_all_set(:,:,t_Nr) = f_true_matrix;

    u_cell = Manipulator_2D_2DoF_get_u_cell(x_all_cell, phi_cell, f_hat_matrix, L1, L2, m1, m2);

    [~, x_all_temp] = ode45( ...
        @(t,x) Manipulator_2D_2DoF_MultiAgent_DynamicFunction(t, x, u_cell, L1, L2, m1, m2), ...
        [t, t+t_step], x_all, opts);
    x_all_next = x_all_temp(end,:)';
    x_all_set(:, t_Nr+1) = x_all_next;
    vartheta_all_set(:, t_Nr+1) = x_all_next - s_all_set(:,t_Nr+1) - ...
        kron(ones(AgentQuantity,1), xl_set(:,t_Nr+1));

    if ~strcmpi(CurrentMode, 'exact')
        %% Online Learning ET: update local GP data for next step
        [LocalGP_set, online_trigger_set, online_trigger_count] = ...
            apply_online_learning_et( ...
            LocalGP_set, online_trigger_set, online_trigger_count, ...
            t_Nr, x_all_matrix, r_matrix, e_cell, ...
            AgentQuantity, y_dim, beta, gamma, ...
            Bidirection_NeighbourSet, eta_underline_set, ...
            MultiAgentSystem, Fl, xi, chi, vartheta_bar);
    end

    % fprintf('t = %6.4f\n', t);
end
TrackingError_vector(end) = norm(vartheta_all_set(:,end));

x_all_matrix_end = reshape(x_all_set(:,end), x_dim, AgentQuantity);
for n = 1:AgentQuantity
    f_true_all_set(:,n,end) = Manipulator_2D_2DoF_UnknownDynamics(x_all_matrix_end(:,n));
    f_hat_all_set(:,n,end)  = f_hat_matrix(:,n);
end

elapsed_time = toc;

fprintf('==================================================');
fprintf('Mode: %s', CurrentMode);
fprintf('Formation: %d', use_formation);
fprintf('Total simulation time: %.2f s', elapsed_time);
fprintf('Final tracking error: %.6f', TrackingError_vector(end));
fprintf('Online learning ET:');
fprintf('  Total triggers: %d', sum(online_trigger_count));
fprintf('  Average triggers: %.2f / agent', mean(online_trigger_count));
fprintf('==================================================');

%% 10. Save
if nargin >= 3
    if ~exist(SaveFolderName,'dir'), mkdir(SaveFolderName); end
    save(fullfile(SaveFolderName,[SaveFileName,'.mat']), ...
        't_set','TrackingError_vector','CurrentMode','use_formation',...
        'f_hat_all_set','f_true_all_set','vartheta_all_set', ...
        'online_trigger_set','online_trigger_count','eta_underline_set', ...
        'vartheta_bar','elapsed_time');
end
end

function P_tp = build_tp_prediction_info( ...
    x_query, LocalGP_set, AgentQuantity, y_dim, method, prior_var, SigmaF)

P_tp = zeros(2*y_dim, AgentQuantity);

for ModelAgentNr = 1:AgentQuantity
    % Use GP of ModelAgentNr, but evaluate at controlled agent's x_query.
    [mu_k, var_k] = LocalGP_set{ModelAgentNr}.predict(x_query);

    mu_k = mu_k(:);
    if isscalar(var_k)
        var_k = var_k * ones(y_dim,1);
    else
        var_k = var_k(:);
    end

    for d = 1:y_dim
        sv = max(var_k(d), 1e-6);
        beta = max(min(0.5*(log(prior_var)-log(sv)), 10), eps);

        switch lower(method)
            case 'moe'
                P_tp(2*d-1,ModelAgentNr) = AgentQuantity * mu_k(d);
                P_tp(2*d,  ModelAgentNr) = AgentQuantity * (sv + mu_k(d)^2);

            case 'gpoe'
                P_tp(2*d-1,ModelAgentNr) = AgentQuantity * beta * mu_k(d) / sv;
                P_tp(2*d,  ModelAgentNr) = AgentQuantity * beta / sv;

            case 'poe'
                P_tp(2*d-1,ModelAgentNr) = AgentQuantity * mu_k(d) / sv;
                P_tp(2*d,  ModelAgentNr) = AgentQuantity / sv;

            case 'bcm'
                P_tp(2*d-1,ModelAgentNr) = AgentQuantity * mu_k(d) / sv;
                P_tp(2*d,  ModelAgentNr) = AgentQuantity / sv - ...
                    (AgentQuantity-1) / prior_var;

            case 'rbcm'
                P_tp(2*d-1,ModelAgentNr) = AgentQuantity * beta * mu_k(d) / sv;
                P_tp(2*d,  ModelAgentNr) = AgentQuantity * beta / sv + ...
                    (1 - AgentQuantity * beta) / prior_var;

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
