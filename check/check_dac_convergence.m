function check_dac_convergence()
rng(0);

%% 参数设置
AgentQuantity = 6; LeaderQuantity = 1;
SystemOrder = 2; q_dim = 2; x_dim = q_dim * SystemOrder;
SigmaF = 1; SigmaL = 0.5*ones(x_dim,1);
SigmaN = 0.05; y_dim = q_dim;
DomainScale = 1.5;
NumInducingPoints = 100;

%% 
AgentConnectionConfigurationSet = {
    {1,2,'Directed',0.5}; {1,6,'Directed',0.5};
    {2,3,'Directed',0.5}; {2,4,'Directed',0.5};
    {3,2,'Directed',0.5}; {3,1,'Directed',0.5};
    {4,3,'Directed',0.5}; {4,5,'Directed',0.5};
    {5,6,'Directed',0.5}; {5,4,'Directed',0.5};
    {6,5,'Directed',0.5}; {6,1,'Directed',0.5}
};
MAS = MultiAgentSystem_Class(6, 1);
MAS.Agent_Topology.addConnectionSet(AgentConnectionConfigurationSet);
MAS.get_ExtendedTopology;
L_lap = MAS.Agent_Topology.LaplacianMatrix;
L_lap = MultiAgentSystem.Agent_Topology.LaplacianMatrix;

%% local GP
LocalGP_set = cell(AgentQuantity,1);
for n = 1:AgentQuantity
    LocalGP_set{n} = LocalGP_MultiOutput(x_dim, y_dim, 100, SigmaN, SigmaF, SigmaL);
    X_in = 2*(rand(x_dim,100)-0.5)*DomainScale;
    Y_in = Manipulator_2D_2DoF_UnknownDynamics(X_in);
    Y_in = Y_in + SigmaN*randn(size(Y_in));
    LocalGP_set{n}.add_Alldata(X_in, Y_in);
    LocalGP_set{n}.xMax = (DomainScale*ones(x_dim,1));
    LocalGP_set{n}.xMin = (-DomainScale*ones(x_dim,1));
end

%% 诱导点
InducingPoints_Coordinates = 2*DomainScale*rand(x_dim,NumInducingPoints) - DomainScale;

%% 计算P（PoE）
[P_inducing, p_dim] = gp_masked_aggregation_init(...
    LocalGP_set, AgentQuantity, NumInducingPoints, ...
    InducingPoints_Coordinates, 'poe');

%% 
Kappa_P_set = [10, 100, 1000];
T_set = [1, 10, 100, 1000];

fprintf('\n%s\n', repmat('=',1,60));
fprintf('%-10s  %-12s  %-15s\n', 'Kappa_P', 'T_integral', 'max agent Xi diff');
fprintf('%s\n', repmat('-',1,60));

for ki = 1:numel(Kappa_P_set)
    kp = Kappa_P_set(ki);
    for ti = 1:numel(T_set)
        T = T_set(ti);
        Zeta_test = zeros(p_dim, AgentQuantity, NumInducingPoints);
        [~, Zeta_conv] = gp_masked_aggregation_update(...
            P_inducing, Zeta_test, L_lap, kp, AgentQuantity, ...
            NumInducingPoints, T, InducingPoints_Coordinates, ...
            SigmaF, SigmaL, x_dim, 'poe', p_dim);
        Xi_dac = P_inducing - Zeta_conv;
        
        % 检查所有agent的Xi是否一致（收敛到同一个值）
        max_diff = 0;
        for n = 2:AgentQuantity
            diff = max(abs(Xi_dac(:,1,:) - Xi_dac(:,n,:)), [], 'all');
            max_diff = max(max_diff, diff);
        end
        fprintf('%-10d  %-12.1f  %-15.2e\n', kp, T, max_diff);
    end
    fprintf('%s\n', repmat('-',1,60));
end
fprintf('\n');
end