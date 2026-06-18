function [TrackingError_vector, t_set] = run_simulation_inducing_point(CurrentMode, SaveFolderName, SaveFileName, use_formation)
% 新增参数 use_formation（默认true）：false时去掉formation offset，保留leader跟踪
if nargin < 4, use_formation = true; end
rng(0);
%% 1. System Parameters
SystemOrder = 2; q_dim = 2; x_dim = q_dim * SystemOrder;
m1 = 1; m2 = 1; L1 = 1; L2 = 1; g = 9.8;
AgentQuantity = 6; LeaderQuantity = 1;

%% 2. Topology
MultiAgentSystem = Manipulator_2D_2DoF_SetMASTopology(AgentQuantity, LeaderQuantity);
L    = MultiAgentSystem.Agent_Topology.LaplacianMatrix;
B    = MultiAgentSystem.Agent_Leader_Topology.ConnectionMatrix(:,1);

%% 3. Controller Parameters
% 照搬 Test.m 第26-66行
c = 10;
lambda_set = [1; 1];
Fl = 0.25;
lambda_n = lambda_set(end);
lambda_vector = lambda_set(1:SystemOrder-1);
Lambda = [zeros(SystemOrder-2,1), eye(SystemOrder-2);
          -lambda_set(1)/lambda_n, -lambda_set(2:end-1)/lambda_n];
% Pe & Qe
Qes = 1 * eye(SystemOrder-1);
Pes = care(Lambda, [], Qes);
Pe  = kron(Pes, eye(AgentQuantity));
Qe  = kron(Qes, eye(AgentQuantity));
% Pr & Qr
q  = (L + diag(B)) \ ones(AgentQuantity,1);
Pr = diag(1./q);
Qr = Pr*(L+diag(B)) + (L+diag(B))'*Pr;
% Pz & Qz
Pz          = blkdiag(Pr, Pe);
eig_Pz      = eig(Pz);
max_eig_Pz  = max(eig_Pz);
min_eig_Pz  = min(eig_Pz);
t_vec = [zeros(SystemOrder-2,1); 1/lambda_n];
Phi = Pr*kron(lambda_vector'*t_vec, eye(AgentQuantity));
Psi = Pr*kron(lambda_vector'*Lambda, eye(AgentQuantity)) + ...
      kron(t_vec'*Pes, eye(AgentQuantity));
Qz  = [c*lambda_n*Qr-2*Phi, -Psi; -Psi', Qe];
eig_Qz    = eig(Qz);
min_eig_Qz = min(eig_Qz);
if all(real(eig(Qz))>0) && all(real(eig(Lambda))<0)
    fprintf('The controller is stable!\n');
else
    error('Controller is not stable!');
end
% xi & chi — 照搬 Test.m 第64-66行
xi  = 2 * lambda_n / min_eig_Qz * norm(Pr * (L + diag(B)));
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
% 照搬 Test.m 第83-126行（GP参数和offline数据）
SigmaF = 1; SigmaL = 0.5*ones(x_dim,1);
GP_tau = 1e-8; GP_delta = 0.1; y_dim = q_dim;
LocalGP_Quantity = AgentQuantity;
DomainScale = 1.5;
X_min = DomainScale * [-1,-1,-1,-1];
X_max = DomainScale * [ 1, 1, 1, 1];
MaxDataQuantity_set = 600*ones(AgentQuantity,1);
SigmaN_set = 0.01*ones(AgentQuantity,1);

dac_methods = {'poe','gpoe','moe','bcm','rbcm'};
ac_methods  = {'poe_ac','gpoe_ac','moe_ac','bcm_ac','rbcm_ac'};

% offline模式给满数据，其余模式给200个初始数据点
if strcmpi(CurrentMode, 'offline')
    OfflineDataQuantity_set = 1 * MaxDataQuantity_set;
else
    OfflineDataQuantity_set = 200 * ones(AgentQuantity, 1);
end

LocalGP_set = cell(LocalGP_Quantity, 1);
for LocalGP_Nr = 1:LocalGP_Quantity
    MaxDataQuantity   = MaxDataQuantity_set(LocalGP_Nr);
    OfflineDataQuantity = OfflineDataQuantity_set(LocalGP_Nr);
    SigmaN            = SigmaN_set(LocalGP_Nr);
    LocalGP_set{LocalGP_Nr} = LocalGP_MultiOutput(x_dim, y_dim, MaxDataQuantity, ...
        SigmaN, SigmaF, SigmaL);
    X_in = 2*(rand(x_dim, OfflineDataQuantity)-0.5)*DomainScale;
    Y_in = Manipulator_2D_2DoF_UnknownDynamics(X_in);
    Y_in = Y_in + SigmaN*randn(size(Y_in));
    LocalGP_set{LocalGP_Nr}.add_Alldata(X_in, Y_in);
    LocalGP_set{LocalGP_Nr}.tau  = GP_tau;
    LocalGP_set{LocalGP_Nr}.delta = GP_delta;
    LocalGP_set{LocalGP_Nr}.xMax = X_max;
    LocalGP_set{LocalGP_Nr}.xMin = X_min;
end

%% 7. Bidirectional Neighbour Set & Sigma_update — 照搬 Test.m 第128-150行
Bidirection_NeighbourSet        = cell(AgentQuantity, 1);
Sigma_update_aggregation_set    = nan(AgentQuantity, 1);
for AgentNr = 1:AgentQuantity
    AgentNeighbourSet = MultiAgentSystem.Agent_Topology.NeighbourSet{AgentNr};
    for NeighbourNr = numel(AgentNeighbourSet):-1:1
        NeighbourAgentNr = AgentNeighbourSet(NeighbourNr);
        if isempty(find(MultiAgentSystem.Agent_Topology.NeighbourSet{NeighbourAgentNr} == AgentNr, 1))
            AgentNeighbourSet(NeighbourNr) = [];
        end
    end
    Bidirection_NeighbourSet{AgentNr} = AgentNeighbourSet;

    Sigma_update_set    = nan(numel(AgentNeighbourSet)+1, 1);
    Sigma_update_set(1) = LocalGP_set{AgentNr}.SigmaN;
    for Bidirection_NeighbourNr = 1:numel(AgentNeighbourSet)
        Bidirection_NeighbourAgentNr = AgentNeighbourSet(Bidirection_NeighbourNr);
        Sigma_update_set(Bidirection_NeighbourNr+1) = ...
            LocalGP_set{Bidirection_NeighbourAgentNr}.SigmaF;
    end
    Sigma_update_aggregation = sqrt(1 / (sum(Sigma_update_set.^(-2)) / numel(Sigma_update_set)));
    Sigma_update_aggregation_set(AgentNr) = Sigma_update_aggregation;
end

%% 8. ET Parameters — 照搬 Test.m 第151-161行
beta = 0; gamma = 0.0005;
for LocalGP_Nr = 1:LocalGP_Quantity
    [~,~,~,beta_i,~,~] = LocalGP_set{LocalGP_Nr}.predict(zeros(x_dim,1));
    beta = max(beta, beta_i);
end
beta = 1.0 * beta;
eta_underline_set = sqrt(beta) * Sigma_update_aggregation_set + gamma;
vartheta_bar = xi * chi * norm((eye(AgentQuantity) - diag(B)) * ones(AgentQuantity,1) * Fl + eta_underline_set);

%% 9. Inducing Points & Mode Init — 照搬原始文件第77-115行
Kappa_P = 1;
NumInducingPoints = 400;
InducingPoints_Coordinates = 2*DomainScale*rand(x_dim, NumInducingPoints) - DomainScale;
L_lap = MultiAgentSystem.Agent_Topology.LaplacianMatrix;

if ismember(lower(CurrentMode), dac_methods)
    base_method = lower(CurrentMode);
    [P_inducing, p_dim] = gp_masked_aggregation_init( ...
        LocalGP_set, AgentQuantity, NumInducingPoints, ...
        InducingPoints_Coordinates, base_method);
    Zeta_vector_inducing  = zeros(p_dim, AgentQuantity, NumInducingPoints);
    neighbor_count_per_agent = sum(L_lap < 0, 2);
    Zeta_last_trigger     = P_inducing;
    total_trigger_count   = zeros(AgentQuantity, 1);

    [MaskedGP, Zeta_vector_inducing, Zeta_last_trigger, ~] = gp_masked_aggregation_update( ...
        P_inducing, Zeta_vector_inducing, L_lap, Kappa_P, AgentQuantity, ...
        NumInducingPoints, 0, InducingPoints_Coordinates, SigmaF, SigmaL, ...
        x_dim, base_method, p_dim, Zeta_last_trigger, neighbor_count_per_agent);

elseif ismember(lower(CurrentMode), ac_methods)
    base_method = strrep(lower(CurrentMode), '_ac', '');
    MaskedGP = gp_masked_aggregation_ac( ...
        LocalGP_set, InducingPoints_Coordinates, SigmaF, SigmaL, ...
        x_dim, AgentQuantity, NumInducingPoints, base_method);
end

%% 10. Initial State — 照搬原始文件第117-132行
rng(42);
x_all = rand(x_dim*AgentQuantity, 1);
x_all_set        = nan(x_dim*AgentQuantity, T);
x_all_set(:,1)   = x_all;
vartheta_all_set = nan(x_dim*AgentQuantity, T);
vartheta_all_set(:,1) = x_all - s_all_set(:,1) - kron(ones(AgentQuantity,1), xl_set(:,1));

f_hat_matrix  = zeros(y_dim, AgentQuantity);
f_true_matrix = zeros(y_dim, AgentQuantity);
TrackingError_vector = zeros(1, T);
f_hat_all_set  = nan(y_dim, AgentQuantity, T);
f_true_all_set = nan(y_dim, AgentQuantity, T);
total_trigger_count = zeros(AgentQuantity, 1);

% 照搬 Test.m 第169-176行
mu_cell              = cell(AgentQuantity, AgentQuantity);
var_matrix           = nan(AgentQuantity, AgentQuantity);
eta_matrix           = nan(AgentQuantity, AgentQuantity);
eta_aggregated_vector = nan(AgentQuantity, 1);
trigger_set          = zeros(AgentQuantity, T);

%% 11. Control Loop
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

    % 控制律 — 照搬 Test.m 第194-195行（取回 r_matrix, e_cell 供ET使用）
    [phi_cell, r_matrix, e_cell] = Manipulator_2D_2DoF_ConsensusLaw( ...
        vartheta_cell, x_tilde_cell, x_l_r, MultiAgentSystem, c, lambda_set, s_r_cell);

    % ------------------------------------------------------------------
    % GP预测部分
    % ------------------------------------------------------------------
    switch lower(CurrentMode)
        case dac_methods
            % MaskedGP预测 — 照搬原始文件
            for n = 1:AgentQuantity
                [mu_hat,~] = MaskedGP{n}.predict(x_all_matrix(:,n));
                mu_hat = max(-30, min(30, mu_hat));
                f_hat_matrix(:,n) = mu_hat;
            end

            % DAC Layer ET update — 照搬原始文件
            [MaskedGP_new, Zeta_vector_inducing, Zeta_last_trigger, step_triggers] = ...
                gp_masked_aggregation_update( ...
                P_inducing, Zeta_vector_inducing, L_lap, Kappa_P, AgentQuantity, ...
                NumInducingPoints, t_step, InducingPoints_Coordinates, ...
                SigmaF, SigmaL, x_dim, base_method, p_dim, ...
                Zeta_last_trigger, neighbor_count_per_agent);
            if any(step_triggers)
                MaskedGP = MaskedGP_new;
            end
            total_trigger_count = total_trigger_count + step_triggers;

            % --------------------------------------------------------------
            % Online Learning ET (Learning Layer)
            % 照搬 Test.m 第204-276行
            % --------------------------------------------------------------
            % Step 1: 计算 eta_aggregated_vector — 照搬 Test.m 第204-229行
            for AgentNr = 1:AgentQuantity
                x_i = x_all_matrix(:, AgentNr);
                [mu_cell{AgentNr,AgentNr}, var_matrix(AgentNr,AgentNr), ...
                 eta_matrix(AgentNr,AgentNr)] = ...
                    Manipulator_2D_2DoF_LocalPrediction(x_i, AgentNr, ...
                    LocalGP_set, beta, gamma, y_dim);
                AgentBidirection_NeighbourSet = Bidirection_NeighbourSet{AgentNr};
                for Bidirection_NeighborNr = 1:numel(AgentBidirection_NeighbourSet)
                    Bidirection_NeighborAgentNr = AgentBidirection_NeighbourSet(Bidirection_NeighborNr);
                    [mu_cell{AgentNr, Bidirection_NeighborAgentNr}, ...
                     var_matrix(AgentNr, Bidirection_NeighborAgentNr), ...
                     eta_matrix(AgentNr, Bidirection_NeighborAgentNr)] = ...
                        Manipulator_2D_2DoF_LocalPrediction(x_i, Bidirection_NeighborAgentNr, ...
                        LocalGP_set, beta, gamma, y_dim);
                end
                var_row_vector_SingleAgent = var_matrix(AgentNr, :);
                mu_row_cell_SingleAgent    = mu_cell(AgentNr, :);
                [~, eta_aggregated_i] = ET_MAS_GP_Leader_GPAggregation_SingleAgent( ...
                    AgentNr, AgentBidirection_NeighbourSet, ...
                    var_row_vector_SingleAgent, mu_row_cell_SingleAgent, beta, gamma);
                eta_aggregated_vector(AgentNr) = eta_aggregated_i;
            end

            % Step 2: ET trigger 判断 — 照搬 Test.m 第231-251行
            for AgentNr = 1:AgentQuantity
                trigger_set(AgentNr, t_Nr) = ...
                    Manipulator_2D_2DoF_DistributedET(AgentNr, ...
                    r_matrix, e_cell, eta_underline_set, eta_aggregated_vector, ...
                    MultiAgentSystem, Fl, xi, chi, vartheta_bar);
            end

            % Step 3: Online learning + P重算 — 照搬 Test.m 第252-276行
            for AgentNr = 1:AgentQuantity
                x_i = x_all_matrix(:, AgentNr);
                AgentBidirection_NeighbourSet = Bidirection_NeighbourSet{AgentNr};
                if trigger_set(AgentNr, t_Nr) == 1
                    y_i = Manipulator_2D_2DoF_UnknownDynamics(x_i) + ...
                          LocalGP_set{AgentNr}.SigmaN * randn(y_dim, 1);
                    if LocalGP_set{AgentNr}.DataQuantity >= LocalGP_set{AgentNr}.MaxDataQuantity
                        LocalGP_set{AgentNr}.downdateParam(1);
                    end
                    LocalGP_set{AgentNr}.addPoint(x_i, y_i);

                    % Update Local Prediction — 照搬 Test.m 第264-265行
                    [mu_cell{AgentNr,AgentNr}, var_matrix(AgentNr,AgentNr), ...
                     eta_matrix(AgentNr,AgentNr)] = ...
                        Manipulator_2D_2DoF_LocalPrediction(x_i, AgentNr, ...
                        LocalGP_set, beta, gamma, y_dim);
                    % Update Aggregation — 照搬 Test.m 第267-274行
                    var_row_vector_SingleAgent = var_matrix(AgentNr, :);
                    mu_row_cell_SingleAgent    = mu_cell(AgentNr, :);
                    [~, eta_aggregated_i] = ET_MAS_GP_Leader_GPAggregation_SingleAgent( ...
                        AgentNr, AgentBidirection_NeighbourSet, ...
                        var_row_vector_SingleAgent, mu_row_cell_SingleAgent, beta, gamma);
                    eta_aggregated_vector(AgentNr) = eta_aggregated_i;

                    % 重算 P_inducing 对应列 — share dataset信息进DAC
                    P_inducing = recompute_P_single_agent(P_inducing, AgentNr, ...
                        LocalGP_set, NumInducingPoints, InducingPoints_Coordinates, ...
                        AgentQuantity, base_method);
                end
            end

        case ac_methods
            for n = 1:AgentQuantity
                [mu_hat,~] = MaskedGP{n}.predict(x_all_matrix(:,n));
                mu_hat = max(-30, min(30, mu_hat));
                f_hat_matrix(:,n) = mu_hat;
            end

        case 'local'
            for n = 1:AgentQuantity
                [mu_n,~] = LocalGP_set{n}.predict(x_all_matrix(:,n));
                mu_n = max(-30, min(30, mu_n));
                f_hat_matrix(:,n) = mu_n;
            end

        case 'exact'
            for n = 1:AgentQuantity
                f_hat_matrix(:,n) = Manipulator_2D_2DoF_UnknownDynamics(x_all_matrix(:,n));
            end
    end

    % 记录每步预测值和真实动力学
    for n = 1:AgentQuantity
        f_true_matrix(:,n) = Manipulator_2D_2DoF_UnknownDynamics(x_all_matrix(:,n));
    end
    f_hat_all_set(:,:,t_Nr)  = f_hat_matrix;
    f_true_all_set(:,:,t_Nr) = f_true_matrix;

    % 控制输入
    u_cell = Manipulator_2D_2DoF_get_u_cell(x_all_cell, phi_cell, f_hat_matrix, L1, L2, m1, m2);

    % 仿真
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

% 补最后一步记录
x_all_matrix_end = reshape(x_all_set(:,end), x_dim, AgentQuantity);
for n = 1:AgentQuantity
    f_true_all_set(:,n,end) = Manipulator_2D_2DoF_UnknownDynamics(x_all_matrix_end(:,n));
    f_hat_all_set(:,n,end)  = f_hat_matrix(:,n);
end

fprintf('Mode: %s  Formation: %d  done, total=%.2fs\n', CurrentMode, use_formation, toc);

if ismember(lower(CurrentMode), dac_methods)
    trigger_count_per_agent_point = mean(total_trigger_count) / NumInducingPoints;
    fprintf('DAC ET 平均触发次数: %.2f / agent / inducing point\n', trigger_count_per_agent_point);
    fprintf('Online Learning ET 触发次数:\n');
    disp(sum(trigger_set, 2));
end

%% 12. Save
if ~exist('trigger_count_per_agent_point','var')
    trigger_count_per_agent_point = 0;
end
if nargin >= 3
    if ~exist(SaveFolderName,'dir'), mkdir(SaveFolderName); end
    save(fullfile(SaveFolderName,[SaveFileName,'.mat']), ...
        't_set', 'TrackingError_vector', 'CurrentMode', 'use_formation', ...
        'f_hat_all_set', 'f_true_all_set', 'vartheta_all_set', ...
        'total_trigger_count', 'trigger_count_per_agent_point', ...
        'trigger_set', 'NumInducingPoints', 'Kappa_P');
end
end