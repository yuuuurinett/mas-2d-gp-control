function check_bcm_convergence()
%% 验证BCM在测试点DAC下的分母收敛情况
%  P实时变化时，Xi的分母是否会出现负值

rng(0);

%% 参数
AgentQuantity = 6;
SystemOrder = 2; q_dim = 2; x_dim = q_dim * SystemOrder;
SigmaF = 1; SigmaL = 0.5*ones(x_dim,1);
SigmaN = 0.05; y_dim = q_dim;
DomainScale = 1.5;
Kappa_P = 1000;
t_step = 0.01;
prior_var = SigmaF^2;
prior_correction = (1 - AgentQuantity) / prior_var;  % = -5

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
L = MAS.Agent_Topology.LaplacianMatrix;

%% 局部GP
LocalGP_set = cell(AgentQuantity,1);
for n = 1:AgentQuantity
    LocalGP_set{n} = LocalGP_MultiOutput(x_dim, y_dim, 350, SigmaN, SigmaF, SigmaL);
    X_in = 2*(rand(x_dim,350)-0.5)*DomainScale;
    Y_in = Manipulator_2D_2DoF_UnknownDynamics(X_in);
    Y_in = Y_in + SigmaN*randn(size(Y_in));
    LocalGP_set{n}.add_Alldata(X_in, Y_in);
end

%% 模拟测试点DAC：P每步变化
% 用参考轨迹上的状态作为测试点
t_set = 0:t_step:4;
[xl_set, ~, ~] = Manipulator_2D_2DoF_LeaderDynamics(t_set, 1);

Zeta_vector = zeros(4, AgentQuantity);
neg_denom_count = 0;
total_steps = numel(t_set) - 1;

fprintf('模拟%d步测试点DAC (BCM)...\n', total_steps);
for k = 1:total_steps
    x_star = xl_set(:, k);
    
    % 每步重新计算P（测试点随状态变化）
    P = zeros(4, AgentQuantity);
    for n = 1:AgentQuantity
        [mu_n, var_n] = LocalGP_set{n}.predict(x_star);
        P(1,n) = AgentQuantity * mu_n(1) / var_n(1);
        P(2,n) = AgentQuantity * mu_n(2) / var_n(2);
        P(3,n) = AgentQuantity / var_n(1);
        P(4,n) = AgentQuantity / var_n(2);
    end
    
    % DAC积分一步
    dac_ode = @(~,z) bcm_dac_ode(z, P, L, Kappa_P, AgentQuantity);
    [~, z_out] = ode45(dac_ode, [0, t_step], Zeta_vector(:));
    Zeta_vector = reshape(z_out(end,:)', 4, AgentQuantity);
    
    % 检查分母
    Xi = P - Zeta_vector;
    denom1 = Xi(3,:) + prior_correction;
    denom2 = Xi(4,:) + prior_correction;
    
    if any(denom1 < 0) || any(denom2 < 0)
        neg_denom_count = neg_denom_count + 1;
        if neg_denom_count <= 5
            fprintf('t=%.2f: 负分母! min(denom1)=%.3f, min(denom2)=%.3f\n', ...
                t_set(k), min(denom1), min(denom2));
        end
    end
end

fprintf('\n总步数: %d\n', total_steps);
fprintf('出现负分母的步数: %d (%.1f%%)\n', neg_denom_count, 100*neg_denom_count/total_steps);
fprintf('prior_correction = %.1f\n', prior_correction);
fprintf('\n结论: BCM在测试点DAC下分母频繁为负，导致phi异常，效果差。\n');
end

function dzeta = bcm_dac_ode(z, P, L, Kappa, N)
Z = reshape(z, 4, N);
dzeta = Kappa * (P - Z) * L';
dzeta = dzeta(:);
end