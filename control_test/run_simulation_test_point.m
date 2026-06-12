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
lambda_n = lambda_set(end);
lambda_vector = lambda_set(1:SystemOrder-1);
Lambda = [zeros(SystemOrder-2,1), eye(SystemOrder-2);
          -lambda_set(1)/lambda_n, -lambda_set(2:end-1)/lambda_n];
Qes = eye(SystemOrder-1);
Pes = care(Lambda, [], Qes);
q   = (L + diag(B)) \ ones(AgentQuantity,1);
Pr  = diag(1./q);
Qr  = Pr*(L+diag(B)) + (L+diag(B))'*Pr;
t_vec = [zeros(SystemOrder-2,1); 1/lambda_n];
Phi = Pr*kron(lambda_vector'*t_vec, eye(AgentQuantity));
Psi = Pr*kron(lambda_vector'*Lambda, eye(AgentQuantity)) + ...
      kron(t_vec'*Pes, eye(AgentQuantity));
Qz  = [c*lambda_n*Qr-2*Phi, -Psi; -Psi', kron(Qes,eye(AgentQuantity))];
if ~(all(real(eig(Qz))>0) && all(real(eig(Lambda))<0))
    error('Controller is not stable!');
end

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
GP_tau = 1e-8; GP_delta = 0.01; y_dim = q_dim;
DomainScale = 1.5;
MaxDataQuantity_set     = 400*ones(AgentQuantity,1);
OfflineDataQuantity_set = MaxDataQuantity_set;
SigmaN_set = 0.05*ones(AgentQuantity,1);
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

%% 7. Setup - DAC Zeta 初始化
dac_methods = {'poe','gpoe','moe','bcm','rbcm'};
ac_methods  = {'poe_ac','gpoe_ac','moe_ac','bcm_ac','rbcm_ac'};
Kappa_P = 100;
L_lap   = MultiAgentSystem.Agent_Topology.LaplacianMatrix;

switch lower(CurrentMode)
    case {'poe','gpoe','moe','bcm','local','exact'}
        Zeta_vector = zeros(4, AgentQuantity);
    case 'rbcm'
        Zeta_vector = zeros(6, AgentQuantity);
    otherwise
        Zeta_vector = zeros(4, AgentQuantity);
end

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

    [phi_cell, ~, ~] = Manipulator_2D_2DoF_ConsensusLaw( ...
        vartheta_cell, x_tilde_cell, x_l_r, MultiAgentSystem, c, lambda_set, s_r_cell);

    AgentState_matrix = x_all_matrix;

    if ismember(lower(CurrentMode), dac_methods)
        switch lower(CurrentMode)
            case 'poe'
                [Phi_Xi, Zeta_vector] = gp_test_poe( ...
                    AgentState_matrix, LocalGP_set, L_lap, Kappa_P, AgentQuantity, Zeta_vector, t_step);
            case 'gpoe'
                [Phi_Xi, Zeta_vector] = gp_test_gpoe( ...
                    AgentState_matrix, LocalGP_set, L_lap, Kappa_P, AgentQuantity, Zeta_vector, t_step);
            case 'moe'
                [Phi_Xi, Zeta_vector] = gp_test_moe( ...
                    AgentState_matrix, LocalGP_set, L_lap, Kappa_P, AgentQuantity, Zeta_vector, t_step);
            case 'bcm'
                [Phi_Xi, Zeta_vector] = gp_test_bcm( ...
                    AgentState_matrix, LocalGP_set, L_lap, Kappa_P, AgentQuantity, Zeta_vector, t_step);
            case 'rbcm'
                [Phi_Xi, Zeta_vector] = gp_test_rbcm( ...
                    AgentState_matrix, LocalGP_set, L_lap, Kappa_P, AgentQuantity, Zeta_vector, t_step);
        end
        f_hat_matrix = Phi_Xi;

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

    fprintf('t = %6.4f\n', t);
end
TrackingError_vector(end) = norm(vartheta_all_set(:,end));

x_all_matrix_end = reshape(x_all_set(:,end), x_dim, AgentQuantity);
for n = 1:AgentQuantity
    f_true_all_set(:,n,end) = Manipulator_2D_2DoF_UnknownDynamics(x_all_matrix_end(:,n));
    f_hat_all_set(:,n,end)  = f_hat_matrix(:,n);
end

fprintf('Mode: %s  Formation: %d  done, total=%.2fs\n', CurrentMode, use_formation, toc);

%% 10. Save
if nargin >= 3
    if ~exist(SaveFolderName,'dir'), mkdir(SaveFolderName); end
    save(fullfile(SaveFolderName,[SaveFileName,'.mat']), ...
        't_set','TrackingError_vector','CurrentMode','use_formation',...
        'f_hat_all_set','f_true_all_set','vartheta_all_set');
end
end