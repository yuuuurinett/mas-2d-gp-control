function [TrackingError_vector, t_set] = run_simulation_cen(CurrentMode, SaveFolderName, SaveFileName, use_formation)
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
t_start = 0; t_end = 4; t_step = 0.01;
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

%% 7. Initial State
rng(42);  % 固定初始状态seed，确保所有方法一致
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

%% 8. Control Loop
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

    for n = 1:AgentQuantity
        x_n = x_all_matrix(:, n);
        mu_all  = zeros(y_dim, AgentQuantity);
        var_all = zeros(y_dim, AgentQuantity);
        for k = 1:AgentQuantity
            [mu_k, var_k] = LocalGP_set{k}.predict(x_n);
            mu_all(:,k)  = mu_k;
            var_all(:,k) = var_k;
        end
        mu_agg = zeros(y_dim, 1);
        for d = 1:y_dim
            mu_a  = mu_all(d,:)';
            var_a = var_all(d,:)';
            beta  = max(0.5*(log(prior_var) - log(var_a)), eps);
            switch lower(CurrentMode)
                case 'moe'
                    mu_agg(d) = mean(mu_a);
                case 'gpoe'
                    prec = sum(beta./var_a);
                    mu_agg(d) = sum(beta.*mu_a./var_a) / max(prec, eps);
                case 'poe'
                    prec = sum(1./var_a);
                    mu_agg(d) = sum(mu_a./var_a) / max(prec, eps);
                case 'bcm'
                    prec = sum(1./var_a) - (AgentQuantity-1)/prior_var;
                    mu_agg(d) = sum(mu_a./var_a) / max(prec, eps);
                case 'rbcm'
                    prec = sum(beta./var_a) + (1-sum(beta))/prior_var;
                    mu_agg(d) = sum(beta.*mu_a./var_a) / max(prec, eps);
            end
        end
        mu_agg = max(-30, min(30, mu_agg));
        f_hat_matrix(:, n) = mu_agg;
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

fprintf('Mode: CEN-%s  Formation: %d  done, total=%.2fs\n', CurrentMode, use_formation, toc);

%% 9. Save
if nargin >= 3
    if ~exist(SaveFolderName,'dir'), mkdir(SaveFolderName); end
    save(fullfile(SaveFolderName,[SaveFileName,'.mat']), ...
        't_set','TrackingError_vector','CurrentMode','use_formation',...
        'f_hat_all_set','f_true_all_set','vartheta_all_set');
end
end