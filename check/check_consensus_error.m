function check_consensus_error()
%% 检查不同拓扑下DAC的consensus error
%  比较：导师拓扑 vs 强连通拓扑
%  方法：PoE（最简单）
%  指标：所有agent的Xi是否收敛到同一个值

rng(0);

%% 参数
AgentQuantity = 6;
SystemOrder = 2; q_dim = 2; x_dim = q_dim * SystemOrder;
SigmaF = 1; SigmaL = 0.5*ones(x_dim,1);
SigmaN = 0.05; y_dim = q_dim;
DomainScale = 1.5;
NumInducingPoints = 100;
t_step = 0.01;
Kappa_P = 1000;

%% 两种拓扑
% 拓扑1：
MAS1 = Manipulator_2D_2DoF_SetMASTopology(AgentQuantity, 1);
L1 = MAS1.Agent_Topology.LaplacianMatrix;
fprintf('导师拓扑零特征值个数: %d\n', sum(abs(real(eig(L1))) < 1e-6));

% 拓扑2：强连通拓扑
AgentConnectionConfigurationSet2 = {
    {1,2,'Directed',0.5}; {1,6,'Directed',0.5};
    {2,3,'Directed',0.5}; {2,4,'Directed',0.5};
    {3,2,'Directed',0.5}; {3,1,'Directed',0.5};
    {4,3,'Directed',0.5}; {4,5,'Directed',0.5};
    {5,6,'Directed',0.5}; {5,4,'Directed',0.5};
    {6,5,'Directed',0.5}; {6,1,'Directed',0.5}
};
MAS2 = MultiAgentSystem_Class(AgentQuantity, 1);
MAS2.Agent_Topology.addConnectionSet(AgentConnectionConfigurationSet2);
MAS2.get_ExtendedTopology;
L2 = MAS2.Agent_Topology.LaplacianMatrix;
fprintf('强连通拓扑零特征值个数: %d\n\n', sum(abs(real(eig(L2))) < 1e-6));

%% 局部GP
LocalGP_set = cell(AgentQuantity,1);
for n = 1:AgentQuantity
    LocalGP_set{n} = LocalGP_MultiOutput(x_dim, y_dim, 100, SigmaN, SigmaF, SigmaL);
    X_in = 2*(rand(x_dim,100)-0.5)*DomainScale;
    Y_in = Manipulator_2D_2DoF_UnknownDynamics(X_in);
    Y_in = Y_in + SigmaN*randn(size(Y_in));
    LocalGP_set{n}.add_Alldata(X_in, Y_in);
end

%% 诱导点和P
InducingPoints_Coordinates = 2*DomainScale*rand(x_dim,NumInducingPoints) - DomainScale;
[P_inducing, p_dim] = gp_masked_aggregation_init(...
    LocalGP_set, AgentQuantity, NumInducingPoints, ...
    InducingPoints_Coordinates, 'poe');

%% AC参考值：mean(P) 直接算全局均值
Xi_ac = repmat(mean(P_inducing, 2), 1, AgentQuantity, 1);  % 收敛目标

%% 模拟仿真步数累积，记录consensus error
Steps_set = [10, 50, 100, 200, 400];

fprintf('%s\n', repmat('=',1,55));
fprintf('%-20s  %-10s  %-15s\n', 'Topology', 'Steps', 'Consensus Error');
fprintf('%s\n', repmat('-',1,55));

for topology_idx = 1:2
    if topology_idx == 1
        L = L1; topo_name = 'PhD Topology';
    else
        L = L2; topo_name = 'Strongly Connected';
    end
    
    Zeta = zeros(p_dim, AgentQuantity, NumInducingPoints);
    step_idx = 1;
    
    for k = 1:Steps_set(end)
        [~, Zeta] = gp_masked_aggregation_update(...
            P_inducing, Zeta, L, Kappa_P, AgentQuantity, ...
            NumInducingPoints, t_step, InducingPoints_Coordinates, ...
            SigmaF, SigmaL, x_dim, 'poe', p_dim);
        
        if k == Steps_set(step_idx)
            Xi = P_inducing - Zeta;
            % Consensus error: 所有agent的Xi与均值的偏差
            Xi_mean = mean(Xi, 2);
            consensus_err = max(abs(Xi - repmat(Xi_mean, 1, AgentQuantity, 1)), [], 'all');
            fprintf('%-20s  %-10d  %-15.2e\n', topo_name, k, consensus_err);
            step_idx = step_idx + 1;
            if step_idx > numel(Steps_set), break; end
        end
    end
    fprintf('%s\n', repmat('-',1,55));
end

%% AC的consensus error（理论上为0，因为直接取均值）
fprintf('%-20s  %-10s  %-15.2e\n', 'AC (direct mean)', 'N/A', 0.0);
fprintf('%s\n', repmat('=',1,55));
fprintf('\n结论：\n');
fprintf('- Consensus error接近0 → DAC收敛，所有agent达成一致\n');
fprintf('- Consensus error很大 → DAC没收敛\n');
end