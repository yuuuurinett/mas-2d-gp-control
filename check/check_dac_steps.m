function check_dac_steps()
%% 验证DAC在仿真步数累积下的收敛性
%  模拟实际仿真中每步积分t_step=0.01的累积效果

rng(0);

%% 参数设置
AgentQuantity = 6;
SystemOrder = 2; q_dim = 2; x_dim = q_dim * SystemOrder;
SigmaF = 1; SigmaL = 0.5*ones(x_dim,1);
SigmaN = 0.05; y_dim = q_dim;
DomainScale = 1.5;
NumInducingPoints = 100;
t_step = 0.01;

%% 强连通拓扑
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

%% 局部GP
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

%% 诱导点和P
InducingPoints_Coordinates = 2*DomainScale*rand(x_dim,NumInducingPoints) - DomainScale;
[P_inducing, p_dim] = gp_masked_aggregation_init(...
    LocalGP_set, AgentQuantity, NumInducingPoints, ...
    InducingPoints_Coordinates, 'poe');

%% 不同Kappa_P，模拟累积步数
Kappa_P_set = [10, 100, 1000];
Steps_set = [10, 50, 100, 200, 400];  % 对应仿真时间 0.1, 0.5, 1, 2, 4秒

fprintf('\n%s\n', repmat('=',1,55));
fprintf('%-10s  %-10s  %-15s\n', 'Kappa_P', 'Steps', 'max agent Xi diff');
fprintf('%s\n', repmat('-',1,55));

for ki = 1:numel(Kappa_P_set)
    kp = Kappa_P_set(ki);
    Zeta = zeros(p_dim, AgentQuantity, NumInducingPoints);
    step_idx = 1;
    
    for k = 1:Steps_set(end)
        [~, Zeta] = gp_masked_aggregation_update(...
            P_inducing, Zeta, L_lap, kp, AgentQuantity, ...
            NumInducingPoints, t_step, InducingPoints_Coordinates, ...
            SigmaF, SigmaL, x_dim, 'poe', p_dim);
        
        if k == Steps_set(step_idx)
            Xi = P_inducing - Zeta;
            max_diff = 0;
            for n = 2:AgentQuantity
                diff = max(abs(Xi(:,1,:) - Xi(:,n,:)), [], 'all');
                max_diff = max(max_diff, diff);
            end
            fprintf('%-10d  %-10d  %-15.2e\n', kp, k, max_diff);
            step_idx = step_idx + 1;
            if step_idx > numel(Steps_set), break; end
        end
    end
    fprintf('%s\n', repmat('-',1,55));
end
fprintf('\n');
end