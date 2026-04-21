clc;
clear; close all;
rng(0);
%% Multi Agent System Topology
AgentQuantity=6;
AgentConnectionConfigurationSet = {
    
    {1,2,'Directed',0.5};
    {1,6,'Directed',0.5};
    {2,3,'Directed',0.5};
    {2,4,'Directed',0.5};  
    {3,2,'Directed',0.5};
    {3,1,'Directed',0.5};
    {4,3,'Directed',0.5};
    {4,5,'Directed',0.5};
    {5,6,'Directed',0.5};
    {5,4,'Directed',0.5};
    {6,5,'Directed',0.5};
    {6,1,'Directed',0.5};

    };
   
%{
    {1,2,'Directed',1};
    {1,6,'Directed',1};
    {2,3,'Directed',1};
    {3,2,'Directed',1};
    {4,3,'Directed',1};
    {4,5,'Directed',1};
    {5,6,'Directed',1};
    {6,5,'Directed',1};
    {3,1,'Directed',1};  
    {5,4,'Directed',1}; 
%}
 
%% Agent Leader Topology
LeaderQuantity = 1;
AgentLeaderConnectionConfigurationSet = {
    {'Leader',1,'Agent',1,1};
    % {'Leader',1,'Agent',2,1};
    % {'Leader',1,'Agent',3,1};
    {'Leader',1,'Agent',4,1}
    };
%% Create Multi-Agent System Object
MultiAgentSystem = MultiAgentSystem_Class(AgentQuantity,LeaderQuantity);
%% Configure Topology
% Agent Topology
MultiAgentSystem.Agent_Topology.addConnectionSet(AgentConnectionConfigurationSet);
MultiAgentSystem.Agent_Topology.check_Connectivity(true);
% Agent Leader Tolology
MultiAgentSystem.Agent_Leader_Topology.addConnectionSet( ...
    AgentLeaderConnectionConfigurationSet);
%
MultiAgentSystem.get_ExtendedTopology;
%
L = MultiAgentSystem.Agent_Topology.LaplacianMatrix;
B = MultiAgentSystem.Agent_Leader_Topology.ConnectionMatrix;

%% Time Horizon
t_start = 0;
t_end = 4;
t_step = 0.01;
t_set = t_start:t_step:t_end;
%% Set System Dimension
SystemOrder = 2;
q_dim = 2;
x_dim = q_dim * SystemOrder;
%% Set Model Paramter
m1 = 1; m2 = 1;
L1 = 1; L2 = 1;
g  = 9.8;
%% Set Reference Trajectory
x_trajectory_set = nan(x_dim, numel(t_set));
x_trajectory_set(1, :) = sin(0.5*t_set);      % x_{l,1,1}
x_trajectory_set(2, :) = cos(0.5*t_set);      % x_{l,1,2}
x_trajectory_set(3, :) = 0.5*cos(0.5*t_set);  % x_{l,2,1}
x_trajectory_set(4, :) = -0.5*sin(0.5*t_set); % x_{l,2,2}

%% Set Gaussian Processes
% Common Setting
SigmaF = 1;
SigmaL = 0.5 * ones(x_dim,1);
GP_tau = 1e-8;
GP_delta = 0.1;
y_dim = q_dim;
LocalGP_Quantity = AgentQuantity;
DomainScale = 1.5;
X_min = DomainScale * [-1, -1, -1, -1];
X_max = DomainScale * [ 1,  1,  1,  1];

% D1(左下), D2(左上), D3(右上), D4(右下)
Local_X_min = DomainScale * [-1, -1, -1, -1;
    -1, -1,  0,  0;
    0,  0,  0,  0;
    0,  0, -1, -1;
    -1,  0, -1,  0;
    0, -1,  0, -1];

Local_X_max = DomainScale * [ 0,  0,  0,  0;
    0,  0,  1,  1;
    1,  1,  1,  1;
    1,  1,  0,  0;
    0,  1,  0,  1;
    1,  0,  1,  0];

MaxDataQuantity_set = 100 * ones(AgentQuantity,1);
OfflineDataQuantity_set = 100 * ones(AgentQuantity,1);
SigmaN_set = 0.1 * ones(AgentQuantity,1);

% Set Local GP
LocalGP_set = cell(LocalGP_Quantity,1);
for LocalGP_Nr = 1:LocalGP_Quantity
    MaxDataQuantity = MaxDataQuantity_set(LocalGP_Nr);
    OfflineDataQuantity = OfflineDataQuantity_set(LocalGP_Nr);
    SigmaN = SigmaN_set(LocalGP_Nr);

    LocalGP_set{LocalGP_Nr} = LocalGP_MultiOutput(x_dim,y_dim, ...
        MaxDataQuantity,SigmaN,SigmaF,SigmaL);
  
    Data_mu = (Local_X_max(LocalGP_Nr,:) + Local_X_min(LocalGP_Nr,:)) / 2;
    Data_Le = (Local_X_max(LocalGP_Nr,:) - Local_X_min(LocalGP_Nr,:)) / 2;
    X_in = 2 * (rand(x_dim,OfflineDataQuantity) - 0.5) .* ... 
    repmat(Data_Le',[1,OfflineDataQuantity]) + Data_mu';

    Y_in = zeros(y_dim, OfflineDataQuantity);
    for d = 1:OfflineDataQuantity
        Y_in(:,d) = Manipulator_2D_2DoF_UnknownDynamics(X_in(:,d));
    end
    Y_in = Y_in + SigmaN * randn(size(Y_in));

    LocalGP_set{LocalGP_Nr}.add_Alldata(X_in, Y_in); 
    LocalGP_set{LocalGP_Nr}.tau = GP_tau;
    LocalGP_set{LocalGP_Nr}.delta = GP_delta;
    LocalGP_set{LocalGP_Nr}.xMax = max(X_in, [], 2);
    LocalGP_set{LocalGP_Nr}.xMin = min(X_in, [], 2);
end

%% simulation setup
rng(0);
AgentInitialState_matrix = rand(4,AgentQuantity) - 0.5;
lambda_set = [7/4; 1];
Kappa_C = 100; %control gain 

%SimulationModes = {'distributed_POE', 'masked_POE', 'distributed_RBCM','local', 'none', 'exact'};
SimulationModes = {'masked_POE'};
NumModes = numel(SimulationModes);
TrackingError_matrix = zeros(NumModes, numel(t_set));
AdjMatrix  = MultiAgentSystem.Agent_Topology.AdjacencyMatrix;
NeighbourSet = MultiAgentSystem.Agent_Topology.NeighbourSet;

%inducing points initialization
NumInducingPoints = 100;
InducingPoints_Coordinates = 2 * DomainScale * ...
    rand(x_dim, NumInducingPoints) - DomainScale;

P = gp_masked_poe_init(LocalGP_set, AgentQuantity, ...
    NumInducingPoints, InducingPoints_Coordinates);

MaskedGP_Y_diff_history = zeros(1, numel(t_set));

%% control loop
for ModeNr = 1:NumModes
    CurrentMode = SimulationModes{ModeNr}; 
    tic;
    t_gp = 0; t_ode = 0;
    AgentState_matrix = AgentInitialState_matrix;
  

    switch CurrentMode
        case 'distributed_POE'
            Kappa_P = 1000;  
            Zeta_vector = zeros(4, AgentQuantity);

        case 'masked_POE'
            Kappa_P = 10;
            Zeta_vector = zeros(4, AgentQuantity); 
            
            % 初始化诱导点和 MaskedGP
            P = gp_masked_poe_init(LocalGP_set, AgentQuantity, NumInducingPoints, InducingPoints_Coordinates);
            Zeta_vector_inducing = zeros(4, AgentQuantity, NumInducingPoints);
            [MaskedGP_POE, Zeta_vector_inducing] = gp_masked_poe_update( ...
                P, Zeta_vector_inducing, L, Kappa_P, AgentQuantity, ...
                NumInducingPoints, 0, InducingPoints_Coordinates, SigmaF, SigmaL, x_dim);      
            
        case 'distributed_RBCM'
            Kappa_P = 1000;  
            Zeta_vector = zeros(6, AgentQuantity); 
            
        otherwise % 包含 'local', 'exact', 'none'
            Kappa_P = 0;  
            Zeta_vector = zeros(4, AgentQuantity);
    end

   
    for k = 1:numel(t_set)
        CurrentTime = t_set(k);
        LeaderState_vector = x_trajectory_set(:, k);
        TrackingError_matrix(ModeNr, k) = norm(AgentState_matrix(:) ...
            - repmat(LeaderState_vector, AgentQuantity, 1));
        if k == numel(t_set), break; end

        Phi_Xi_vector = zeros(q_dim, AgentQuantity);
        tic_gp = tic;

        switch CurrentMode
            case 'distributed_POE'
                [Phi_Xi_vector, Zeta_vector, Xi_diff] = gp_poe( ...
                    AgentState_matrix, LocalGP_set, L,  ...
                    Kappa_P, AgentQuantity, Zeta_vector, t_step);

                if mod(k, 50) == 0
                    fprintf('k=%d, t=%.2f, Xi diff = %.4e\n', k, t_set(k), Xi_diff);
                end
  
            case 'masked_POE'
                for AgentNr = 1:AgentQuantity
                    x_i = AgentState_matrix(:, AgentNr);
                    [mu_hat, ~] = MaskedGP_POE{AgentNr}.predict(x_i);
                    Phi_Xi_vector(:, AgentNr) = mu_hat;
                end
                % Update DAC and rebuild MaskedGP
                [MaskedGP_POE, Zeta_vector_inducing] = gp_masked_poe_update( ...
                    P, Zeta_vector_inducing,  L, ...
                    Kappa_P, AgentQuantity, NumInducingPoints, t_step, ...
                    InducingPoints_Coordinates, SigmaF, SigmaL, x_dim);

                % 记录收敛情况
                Y_diff_total = 0;
                for AgentNr_i = 1:AgentQuantity
                    for AgentNr_j = AgentNr_i+1:AgentQuantity
                        Y_diff_total = Y_diff_total + ...
                            norm(MaskedGP_POE{AgentNr_i}.Y - MaskedGP_POE{AgentNr_j}.Y, 'fro');
                    end
                end
                MaskedGP_Y_diff_history(k) = Y_diff_total;


            case 'distributed_RBCM'
                [Phi_Xi_vector, Zeta_vector] = gp_rbcm( ...
                    AgentState_matrix, LocalGP_set, ...                  
                    Kappa_P, AgentQuantity, Zeta_vector, t_step);

            case 'local'
                for AgentNr = 1:AgentQuantity
                    [mu_n, ~] = LocalGP_set{AgentNr}.predict(AgentState_matrix(:, AgentNr));
                    Phi_Xi_vector(:, AgentNr) = mu_n;
                end

            case 'none'
                Phi_Xi_vector = zeros(q_dim, AgentQuantity);

            case 'exact'
                for AgentNr = 1:AgentQuantity
                    Phi_Xi_vector(:, AgentNr) = Manipulator_2D_2DoF_UnknownDynamics( ...
                        AgentState_matrix(:, AgentNr));
                end
        end


        t_gp = t_gp + toc(tic_gp);

        MultiAgent_Dynamics_Handle = ...
            @(t, x) multi_agent_dynamics( ...
            t, x, Phi_Xi_vector, AdjMatrix, NeighbourSet, B, lambda_set, ...
            Kappa_C, AgentQuantity, SystemOrder, q_dim, L1, L2, m1, m2);

        tic_ode = tic;
        opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
        [~, x_output] = ode45(MultiAgent_Dynamics_Handle, ...
            [CurrentTime, CurrentTime + t_step], AgentState_matrix(:), opts);
        AgentState_matrix = reshape(x_output(end,:)', x_dim, AgentQuantity);
        t_ode = t_ode + toc(tic_ode);

        fprintf('t = %6.4f\n', CurrentTime);
    end
    fprintf('Mode: %s done, 总时间=%.2f s, GP=%.2f s, ODE=%.2f s\n', ...
        CurrentMode, toc, t_gp, t_ode);

    if strcmp(CurrentMode, 'masked_POE')
    figure('Color','w');
    semilogy(t_set(1:end-1), MaskedGP_Y_diff_history(1:end-1), ...
        'g-', 'LineWidth', 1.5);
    xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\sum_{i<j} \|Y_i - Y_j\|_F$', 'Interpreter', 'latex', 'FontSize', 12);
    title('MaskedGP Convergence across All Agent Pairs');
    xlim([t_start, t_end]);
    grid on; box on;
    set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
    end

end


%% compute ultimate error bound
L_tilde = MultiAgentSystem.Extended_Topology.LaplacianMatrix ...
(1:AgentQuantity, 1:AgentQuantity);

sigma_L_tilde     = min(svd(L_tilde));
sigma_bar_L_tilde = max(svd(L_tilde));

LyapunovQ_gain = 600;
Q_val          = LyapunovQ_gain;
P_val          = Q_val / (2 * lambda_set(1));
Lambda_norm    = lambda_set(1);

% (式46,47)
iota = norm(eye(AgentQuantity) - L_tilde) * lambda_set(1);
Upsilon = [Kappa_C * sigma_L_tilde - iota, ...
    -0.5*(1 + lambda_set(1)) - 0.5*sqrt(P_val);
    -0.5*(1 + lambda_set(1)) - 0.5*sqrt(P_val), ...
    0.5*Q_val];
sigma_Upsilon = min(eig(Upsilon));

% (式49)
C = sqrt(AgentQuantity) * sigma_bar_L_tilde * sqrt(P_val) * (1 + Lambda_norm) / ...
    (sigma_Upsilon * sigma_L_tilde);

% Delta_hat
Delta_hat = 1;

% beta (式14)
GP_tau_bound = 0.01;
r_Omega  = 2*sqrt(2);
beta_val = 2*2*log(r_Omega*sqrt(2)) - log(2*GP_tau_bound) - 2*log(GP_delta/AgentQuantity);
beta_val = max(beta_val, 0);

% 沿 leader 轨迹计算 v(t)
bound_distributed  = zeros(1, numel(t_set));
bound_local = zeros(1, numel(t_set));
bound_exact = zeros(1, numel(t_set));

for k = 1:numel(t_set)
    LeaderState_k = x_trajectory_set(:, k);

    sigma_sum_dist  = 0;
    sigma_sum_local = 0;
    precison_sum        = 0;

    for AgentNr = 1:AgentQuantity
        [~, var_n] = LocalGP_set{AgentNr}.predict(LeaderState_k);
        var_n = var_n(1);
        sigma_sum_dist  = sigma_sum_dist  + (1/var_n) * sqrt(var_n);
        sigma_sum_local = sigma_sum_local + sqrt(var_n);
        precison_sum        = precison_sum + 1/var_n;
    end

    eta_GP_dist  = 2*sqrt(beta_val) * sigma_sum_dist / precison_sum;
    eta_GP_local = 2*sqrt(beta_val) * sigma_sum_local / AgentQuantity;

    bound_distributed(k)  = C * (eta_GP_dist  + Delta_hat);
    bound_local(k) = C * (eta_GP_local + Delta_hat);
    bound_exact(k) = C * Delta_hat;
end


%% Plot
figure('Color','w','Position',[100 100 700 550]);

% 上图: ||e||
subplot(2,1,1);
styles = {'r-','g-','m-','k-','b-','b--'};
hold on; grid on; box on;
set(gca,'YScale','log','FontSize',11,'FontName','Times New Roman');
for ModeNr = 1:NumModes
    plot(t_set, TrackingError_matrix(ModeNr,:), styles{ModeNr}, 'LineWidth',1.5);
end
ylabel('$\|e\|$','Interpreter','latex','FontSize',12);
xlim([t_start, t_end]);
legend('distributed POE','masked POE','distributed RBCM','local','none','exact','Location','northeast','FontSize',10);
set(gca,'XTickLabel',[]);

% 下图: v(t)
subplot(2,1,2);
hold on; grid on; box on;
set(gca,'FontSize',11,'FontName','Times New Roman');
plot(t_set, bound_local, 'k-',  'LineWidth',1.5, 'DisplayName','local');
plot(t_set, bound_distributed,  'r-',  'LineWidth',1.5, 'DisplayName','distributed');
plot(t_set, bound_exact, 'b--', 'LineWidth',1.5, 'DisplayName','exact');
ylabel('$v(t)$','Interpreter','latex','FontSize',12);
xlabel('$t$','Interpreter','latex','FontSize',12);
xlim([t_start, t_end]);
legend('Location','northeast','FontSize',10);
