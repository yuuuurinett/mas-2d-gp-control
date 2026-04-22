rng(0); 
%% 1. Multi Agent System Topology
AgentQuantity = 6;
AgentConnectionConfigurationSet = {
    {1,2,'Directed',0.5}; {1,6,'Directed',0.5};
    {2,3,'Directed',0.5}; {2,4,'Directed',0.5};  
    {3,2,'Directed',0.5}; {3,1,'Directed',0.5};
    {4,3,'Directed',0.5}; {4,5,'Directed',0.5};
    {5,6,'Directed',0.5}; {5,4,'Directed',0.5};
    {6,5,'Directed',0.5}; {6,1,'Directed',0.5}
};

LeaderQuantity = 1;
AgentLeaderConnectionConfigurationSet = {
    {'Leader',1,'Agent',1,1};
    {'Leader',1,'Agent',4,1}
};

MultiAgentSystem = MultiAgentSystem_Class(AgentQuantity,LeaderQuantity);
MultiAgentSystem.Agent_Topology.addConnectionSet(AgentConnectionConfigurationSet);
MultiAgentSystem.Agent_Topology.check_Connectivity(true);
MultiAgentSystem.Agent_Leader_Topology.addConnectionSet(AgentLeaderConnectionConfigurationSet);
MultiAgentSystem.get_ExtendedTopology;

L = MultiAgentSystem.Agent_Topology.LaplacianMatrix;
B = MultiAgentSystem.Agent_Leader_Topology.ConnectionMatrix;

%% 2. Time Horizon & Parameters
t_start = 0; t_end = 4; t_step = 0.01;
t_set = t_start:t_step:t_end;

SystemOrder = 2; q_dim = 2; x_dim = q_dim * SystemOrder;
m1 = 1; m2 = 1; L1 = 1; L2 = 1; g = 9.8;

%% 3. Set Reference Trajectory
x_trajectory_set = nan(x_dim, numel(t_set));
x_trajectory_set(1, :) = sin(0.5*t_set);
x_trajectory_set(2, :) = cos(0.5*t_set);
x_trajectory_set(3, :) = 0.5*cos(0.5*t_set);
x_trajectory_set(4, :) = -0.5*sin(0.5*t_set);

%% 4. Set Gaussian Processes
SigmaF = 1; SigmaL = 0.5 * ones(x_dim,1);
GP_tau = 1e-8; GP_delta = 0.1; y_dim = q_dim;
LocalGP_Quantity = AgentQuantity;
DomainScale = 1.5;

Local_X_min = DomainScale * [-1, -1, -1, -1; -1, -1,  0,  0; 0,  0,  0,  0; 0,  0, -1, -1; -1,  0, -1,  0; 0, -1,  0, -1];
Local_X_max = DomainScale * [ 0,  0,  0,  0; 0,  0,  1,  1; 1,  1,  1,  1; 1,  1,  0,  0; 0,  1,  0,  1; 1,  0,  1,  0];

MaxDataQuantity_set     = 100 * ones(AgentQuantity,1);
OfflineDataQuantity_set = 100 * ones(AgentQuantity,1);
SigmaN_set              = 0.1 * ones(AgentQuantity,1);

LocalGP_set = cell(LocalGP_Quantity,1);
for LocalGP_Nr = 1:LocalGP_Quantity
    MaxDataQuantity    = MaxDataQuantity_set(LocalGP_Nr);
    OfflineDataQuantity = OfflineDataQuantity_set(LocalGP_Nr);
    SigmaN             = SigmaN_set(LocalGP_Nr);

    LocalGP_set{LocalGP_Nr} = LocalGP_MultiOutput(x_dim, y_dim, MaxDataQuantity, SigmaN, SigmaF, SigmaL);

    Data_mu = (Local_X_max(LocalGP_Nr,:) + Local_X_min(LocalGP_Nr,:)) / 2;
    Data_Le = (Local_X_max(LocalGP_Nr,:) - Local_X_min(LocalGP_Nr,:)) / 2;
    X_in = 2*(rand(x_dim,OfflineDataQuantity)-0.5) .* repmat(Data_Le',[1,OfflineDataQuantity]) + Data_mu';

    Y_in = zeros(y_dim, OfflineDataQuantity);
    for d = 1:OfflineDataQuantity
        Y_in(:,d) = Manipulator_2D_2DoF_UnknownDynamics(X_in(:,d));
    end
    Y_in = Y_in + SigmaN * randn(size(Y_in));

    LocalGP_set{LocalGP_Nr}.add_Alldata(X_in, Y_in);
    LocalGP_set{LocalGP_Nr}.tau   = GP_tau;
    LocalGP_set{LocalGP_Nr}.delta = GP_delta;
    LocalGP_set{LocalGP_Nr}.xMax  = max(X_in, [], 2);
    LocalGP_set{LocalGP_Nr}.xMin  = min(X_in, [], 2);
end

%% 5. Simulation Setup
rng(0);
AgentInitialState_matrix = rand(4, AgentQuantity) - 0.5;
AdjMatrix = MultiAgentSystem.Agent_Topology.AdjacencyMatrix;

%% 6. Test gp_poe_ac
AgentState_matrix = AgentInitialState_matrix;  % <-- 关键：初始化状态

fprintf('=== 验证 gp_poe_ac ===\n');
for k = 1:3
    fprintf('\n--- Step %d ---\n', k);

    % 调用 gp_poe_ac
    Phi_AC = gp_poe_ac(AgentState_matrix, LocalGP_set, AgentQuantity);

    % 验证1: 所有agent的phi是否完全一样
    max_diff = max(max(abs(Phi_AC - Phi_AC(:,1))));
    fprintf('Agent间最大差异: %.2e (应为0)\n', max_diff);

    % 验证2: 和ground truth及local GP比
    for AgentNr = 1:AgentQuantity
        gt = Manipulator_2D_2DoF_UnknownDynamics(AgentState_matrix(:, AgentNr));
        [mu_local, ~] = LocalGP_set{AgentNr}.predict(AgentState_matrix(:, AgentNr));
        fprintf('Agent%d  gt=[%6.3f,%6.3f]  ac=[%6.3f,%6.3f]  local=[%6.3f,%6.3f]\n', ...
            AgentNr, gt(1), gt(2), ...
            Phi_AC(1,AgentNr), Phi_AC(2,AgentNr), ...
            mu_local(1), mu_local(2));
    end
end