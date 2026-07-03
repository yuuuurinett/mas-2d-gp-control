%% ET_MAS_GP_Leader_Test
% clc; 
clear; close all;
% rng(0);
%%
% EventTriggerType = 'distributed'; % centralized, NoTrigger, continuous
EventTriggerType = 'continuous';
%% Multi Agent System Topology
AgentQuantity = 6;
AgentConnectionConfigurationSet = {
    {1, 2, 'Directed', 1}; 
    {1, 6, 'Directed', 1}; 
    {2, 3, 'Undirected', 1};   
    {4, 3, 'Directed', 1}; 
    {4, 5, 'Directed', 1}; 
    {5, 6, 'Undirected', 1}; 
    };
%% Agent Leader Topology
LeaderQuantity = 1;
AgentLeaderConnectionConfigurationSet = {
	{'Leader',1,'Agent',1,1};
% 	{'Leader',1,'Agent',2,1};
% 	{'Leader',1,'Agent',3,1};
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
B = diag(MultiAgentSystem.Agent_Leader_Topology.ConnectionMatrix);
%% Show: Topology
MultiAgentSystem.Agent_Topology.show_topology;
%% Time Horizon
t_start = 0;
t_end = 60;
t_step = 0.02;
t_set = t_start:t_step:t_end;
%% Set System Dimension
SystemOrder = 2;
q_dim = 2;
x_dim = q_dim * SystemOrder;
%% Set Reference Trajectory
[x_l_set, x_l_r_set] = ET_MAS_GP_Leader_LeaderDynamics(t_set);
Fl = 0.01;

s_set_cell = cell(AgentQuantity,1);
for AgentNr = 1:AgentQuantity
	[s_i_set, s_i_r_set] = ET_MAS_GP_Leader_RelativePositionDynamics(t_set, AgentNr);
	s_set_cell{AgentNr}.s_i_set = s_i_set;
	s_set_cell{AgentNr}.s_i_r_set = s_i_r_set;
end

ET_MAS_GP_Leader_Plot_Reference;
%% Control Gain
% lambda
lambda_set = [1;1];
lambda_n = lambda_set(end);
lambda_vector = lambda_set(1:SystemOrder - 1);
Lambda = [	zeros(SystemOrder - 2,1),	eye(SystemOrder - 2);
			- lambda_set(1) / lambda_n,	- lambda_set(2:end-1) / lambda_n];
% Pe & Qe
Qes = 1 * eye(SystemOrder - 1);
Pes = care(Lambda,[],Qes);
Pe = kron(Pes,eye(AgentQuantity));
Qe = kron(Qes,eye(AgentQuantity));
eig_Pe = eig(Pe);
% Pr & Qr
q = (L + B) \ ones(AgentQuantity,1);
Pr = diag(1 ./ q);
Qr = Pr * (L + B) + (L + B)' * Pr;
eig_Pr = eig(Pr);
% Pz & Qz
c = 10;
Pz = blkdiag(Pr,Pe);
eig_Pz = eig(Pz);
max_eig_Pz = max(eig_Pz);
min_eig_Pz = min(eig_Pz);
t = [zeros(SystemOrder - 2,1);1 / lambda_n];
Phi = Pr * kron(lambda_vector' * t, eye(AgentQuantity));
Psi = Pr * kron(lambda_vector' * Lambda, eye(AgentQuantity)) + ...
	kron(t' * Pes,eye(AgentQuantity));
Qz = [c * lambda_n * Qr - 2 * Phi, - Psi; ...
	- Psi', Qe];
eig_Qz = eig(Qz);
min_eig_Qz = min(eig_Qz);
if all(real(eig(Qz)) > 0) && all(real(eig(Lambda)) < 0)
	fprintf('The controller is stable!\n');
else
	error('The controller is not stable!');
end
% xi & chi
xi = 2 * lambda_n / min_eig_Qz * norm(Pr * (L + B));
chi = sqrt((1 + norm([t,Lambda])^2) * max_eig_Pz / min_eig_Pz ) * ...
	norm(inv(L + B));
% chi = sqrt((1 + norm([t,Lambda])^2) ) * norm(inv(L + B));
%% Set Gaussian Processes
% Common Setting
SigmaF = 1;
SigmaL = 0.5 * ones(x_dim,1);
GP_tau = 1e-8;
GP_delta = 0.1;
y_dim = q_dim;
LocalGP_Quantity = AgentQuantity;

DomainScale = 1.5;
X_min = DomainScale * [-1,-1];
X_max = DomainScale * [ 1, 1];
%Local_X_min = DomainScale * [ 0, 0; -1, 0; -1, -1; 0,-1]; % Local_X_min = DomainScale * [-1,-1; -1,-1; -1, -1;-1,-1];
%Local_X_max = DomainScale * [ 1, 1;  0, 1;  0,  0; 1, 0]; % Local_X_max = DomainScale * [ 1, 1;  1, 1;  1,  1; 1, 1];
Local_X_min = repmat(DomainScale * [-1, -1], AgentQuantity, 1);
Local_X_max = repmat(DomainScale * [1, 1], AgentQuantity, 1);

% Individual Setting
MaxDataQuantity_set = 200 * ones(AgentQuantity,1);
%SigmaN_set = [0.01; 0.01; 0.01; 0.01];
SigmaN_set = 0.01 * ones(AgentQuantity, 1);
if strcmpi(EventTriggerType,'distributed') || strcmpi(EventTriggerType,'centralized')
	OfflineDataQuantity_set = 0*[50;50;50;50];
else
	OfflineDataQuantity_set = 4*[50;50;50;50];
end

% Set Local GP
LocalGP_set = cell(LocalGP_Quantity,1);
for LocalGP_Nr = 1:LocalGP_Quantity
	MaxDataQuantity = MaxDataQuantity_set(LocalGP_Nr);
	OfflineDataQuantity = OfflineDataQuantity_set(LocalGP_Nr);
	SigmaN = SigmaN_set(LocalGP_Nr);
	LocalGP_set{LocalGP_Nr} = LocalGP_MultiOutput(x_dim,y_dim,MaxDataQuantity, ...
		SigmaN,SigmaF,SigmaL);

	Data_mu = (Local_X_max(LocalGP_Nr,:) + Local_X_min(LocalGP_Nr,:)) / 2;
	Data_Le = (Local_X_max(LocalGP_Nr,:) - Local_X_min(LocalGP_Nr,:)) / 2;
	X_in = 2 * (rand(x_dim,OfflineDataQuantity) - 0.5) .* repmat(Data_Le',[1,OfflineDataQuantity]) + Data_mu';

	Y_in = ET_MAS_GP_Leader_UnknownDynamics(X_in,q_dim,SystemOrder);
	Y_in = Y_in + SigmaN * randn(size(Y_in));
	LocalGP_set{LocalGP_Nr}.add_Alldata(X_in,Y_in);

	LocalGP_set{LocalGP_Nr}.tau = GP_tau;
	LocalGP_set{LocalGP_Nr}.delta = GP_delta;
	LocalGP_set{LocalGP_Nr}.xMax = X_max;
	LocalGP_set{LocalGP_Nr}.xMin = X_min;
end
%% Plot Initial Offline Training Data Set
ET_MAS_GP_Leader_Plot_GPDataSet;
%%
Bidirection_NeighbourSet = cell(AgentQuantity,1);
Sigma_update_aggregation_set = nan(AgentQuantity,1);
for AgentNr = 1:AgentQuantity
	AgentNeighbourSet = MultiAgentSystem.Agent_Topology.NeighbourSet{AgentNr};
	for NeighbourNr = numel(AgentNeighbourSet):-1:1
		NeighbourAgentNr = AgentNeighbourSet(NeighbourNr);
		if isempty(find(MultiAgentSystem.Agent_Topology.NeighbourSet{NeighbourAgentNr} == AgentNr, 1))
			AgentNeighbourSet(NeighbourNr) = [];
		end
	end
	AgentBidirection_NeighbourSet = AgentNeighbourSet;
	Bidirection_NeighbourSet{AgentNr} = AgentBidirection_NeighbourSet;

	Sigma_update_set = nan(numel(AgentBidirection_NeighbourSet)+1,1);
	Sigma_update_set(1) = LocalGP_set{AgentNr}.SigmaN;
	for Bidirection_NeighbourNr = 1:numel(AgentBidirection_NeighbourSet)
		Bidirection_NeighbourAgentNr = AgentBidirection_NeighbourSet(Bidirection_NeighbourNr);
		Sigma_update_set(Bidirection_NeighbourNr + 1) = LocalGP_set{Bidirection_NeighbourAgentNr}.SigmaF;
	end
	Sigma_update_aggregation = sqrt(1 / (sum(Sigma_update_set .^ (-2)) / numel(Sigma_update_set)));
	Sigma_update_aggregation_set(AgentNr) = Sigma_update_aggregation;
end
% Tracking error bound
beta = 0;
gamma = 0.0005;
for LocalGP_Nr = 1:LocalGP_Quantity
	[~,~,~,beta_i,gamma_i,eta_min] = LocalGP_set{LocalGP_Nr}.predict(zeros(x_dim,1));
% 	[beta_i,gamma_i,eta_min,sqrt(beta_i)*LocalGP_set{LocalGP_Nr}.SigmaN]
	beta = max(beta, beta_i);
end
beta = 0.1 * beta;
eta_underline_set = sqrt(beta) * Sigma_update_aggregation_set + gamma;
vartheta_bar = xi * chi * norm((eye(AgentQuantity) - B) * ones(AgentQuantity,1) * Fl + eta_underline_set);
%% Initial State
x0 = [];
s_set = [];
s_r_set = [];
x_l_0 = x_l_set(:,1);
for AgentNr = 1:AgentQuantity
	s_i_set = s_set_cell{AgentNr}.s_i_set;
	s_i_r_set = s_set_cell{AgentNr}.s_i_r_set;
	s_i_0 = s_i_set(:,1);
	s_set = [s_set;s_i_set];
	s_r_set = [s_r_set;s_i_r_set];
	x0 = [x0; x_l_0 + s_i_0];
end
x0 = 0 * DomainScale * 2 * (rand(AgentQuantity * x_dim,1) - 0.5);
s_0 = s_set(:,1);
x_l_0 = x_l_set(:,1);
%% Simulation
x_set = nan(AgentQuantity * x_dim,numel(t_set));
x_set(:,1) = x0;
if strcmpi(EventTriggerType,'continuous')
	trigger_set = ones(AgentQuantity,numel(t_set));
else
	trigger_set = zeros(AgentQuantity,numel(t_set));
end
vartheta_all_set = nan(AgentQuantity * x_dim,numel(t_set));
vartheta_all_set(:,1) = x0 - s_0 - kron(ones(AgentQuantity,1), x_l_0);
opt = odeset('RelTol',1e-8,'AbsTol',1e-10);
for t_Nr = 1:numel(t_set) - 1
	t = t_set(t_Nr);
	x_all = x_set(:,t_Nr);
	x_matrix = reshape(x_all, [], AgentQuantity);
	s_all = s_set(:,t_Nr);
	s_r_all = s_r_set(:,t_Nr);
	s_r_cell = ET_MAS_GP_Leader_vector2cell(s_r_all, AgentQuantity, 1);
	x_tilde_all = x_all - s_all;
	x_tilde_cell = ET_MAS_GP_Leader_vector2cell(x_tilde_all, AgentQuantity, SystemOrder);
	x_l_r = x_l_r_set(:,t_Nr);
	vartheta_all = vartheta_all_set(:,t_Nr);
	vartheta_cell = ET_MAS_GP_Leader_vector2cell(vartheta_all, AgentQuantity, SystemOrder);

	e_cell = cell(AgentQuantity, SystemOrder);
	u_cell = cell(AgentQuantity,1);
	mu_cell = cell(AgentQuantity, AgentQuantity);
	var_matrix = nan(AgentQuantity, AgentQuantity);
	eta_matrix = nan(AgentQuantity, AgentQuantity);
	eta_aggregated_vector = nan(AgentQuantity,1);
	f_hat_matrix = nan(y_dim, AgentQuantity);
	f_true_matrix = nan(y_dim, AgentQuantity);
	r_matrix = nan(q_dim,AgentQuantity);
	% Consensus Control
	for AgentNr = 1:AgentQuantity
		b_ii = MultiAgentSystem.Agent_Leader_Topology.ConnectionMatrix(AgentNr);
		AgentNeighbourSet = MultiAgentSystem.Agent_Topology.NeighbourSet{AgentNr};
		AgentNeighbourQuantity = numel(AgentNeighbourSet);
		r_i = 0;
		for SystemOrderNr = 1:SystemOrder
			e_ik = b_ii * vartheta_cell{AgentNr, SystemOrderNr};
			x_tilde_ik = x_tilde_cell{AgentNr, SystemOrderNr};
			for NeighborNr = 1:AgentNeighbourQuantity
				NeighborAgentNr = AgentNeighbourSet(NeighborNr);
				a_ij = MultiAgentSystem.Agent_Topology.AdjacencyMatrix(AgentNr, NeighborAgentNr);
				x_tilde_jk = x_tilde_cell{NeighborAgentNr, SystemOrderNr};
				e_ik = e_ik + a_ij * (x_tilde_ik - x_tilde_jk);
			end
			e_cell{AgentNr, SystemOrderNr} = e_ik;
			lambda_k = lambda_set(SystemOrderNr);
			r_i = r_i + lambda_k * e_ik;
		end
		r_matrix(:,AgentNr) = r_i;
		s_r_i = s_r_cell{AgentNr};
		u_cell{AgentNr} = -c * r_i + s_r_i + b_ii * x_l_r;
	end
	% Cooperative Learning
	for AgentNr = 1:AgentQuantity
		b_ii = MultiAgentSystem.Agent_Leader_Topology.ConnectionMatrix(AgentNr);
		% Local prediction
		[mu_cell,var_matrix,eta_matrix] = ET_MAS_GP_Leader_LocalPrediction(AgentNr, ...
			x_matrix,LocalGP_set,beta,gamma,y_dim,mu_cell,var_matrix,eta_matrix);
		% Cooperative Prediction
		AgentBidirection_NeighbourSet = Bidirection_NeighbourSet{AgentNr};
		AgentBidirection_NeighbourQuantity = numel(AgentBidirection_NeighbourSet);
		x_i = x_matrix(:,AgentNr);
		for Bidirection_NeighborNr = 1:AgentBidirection_NeighbourQuantity
			Bidirection_NeighborAgentNr = AgentBidirection_NeighbourSet(Bidirection_NeighborNr);
			[mu_j_xi, var_j_xi] = LocalGP_set{Bidirection_NeighborAgentNr}.predict(x_i);
			eta_j_xi = sqrt(y_dim) * (sqrt(beta * var_j_xi(1)) + gamma);
			mu_cell{AgentNr, Bidirection_NeighborAgentNr} = mu_j_xi;
			var_matrix(AgentNr, Bidirection_NeighborAgentNr) = var_j_xi(1);
			eta_matrix(AgentNr, Bidirection_NeighborAgentNr) = eta_j_xi;
		end
		% GP aggregation - Posterior variance
		var_row_vector_SingleAgent = var_matrix(AgentNr,:);
		mu_row_cell_SingleAgent = mu_cell(AgentNr, :);
		[f_hat_i, eta_aggregated_i] = ET_MAS_GP_Leader_GPAggregation_SingleAgent( ...
			AgentNr, AgentBidirection_NeighbourSet, ...
			var_row_vector_SingleAgent, mu_row_cell_SingleAgent, beta, gamma);
		% Distributed Event-trigger
		if strcmpi(EventTriggerType,'distributed')
			z_i = nan(x_dim,1);
			r_i = r_matrix(:,AgentNr);
			z_i(1:q_dim) = r_i;
			for SystemOrderNr = 1:SystemOrder-1
				z_i(SystemOrderNr * q_dim + (1:q_dim)) = e_cell{AgentNr, SystemOrderNr};
			end
			eta_underline_i = eta_underline_set(AgentNr);
			rho_i = (eta_aggregated_i + (1 - b_ii) * Fl)^2;
			rho_bar_i = xi^(-2) * max(z_i' * z_i - chi^(-2) * vartheta_bar ^ 2 / AgentQuantity, 0) + ...
				(eta_underline_i + (1 - b_ii) * Fl)^2;
			if rho_i > rho_bar_i
				trigger_set(AgentNr, t_Nr) = 1;
			end
		end
		% Save Distributed GP result
		f_hat_matrix(:,AgentNr) = f_hat_i;
		f_true_matrix(:,AgentNr) = MAS_SingleAgent_UnknownDynamics(x_i,q_dim,SystemOrder);
		eta_aggregated_vector(AgentNr) = eta_aggregated_i;
	end
	% Centralized Event-trigger
	if strcmpi(EventTriggerType,'centralized')
		epsilon = nan((SystemOrder - 1) * q_dim * AgentQuantity,1);
		for SystemOrderNr = 1:SystemOrder - 1
			epsilon_k = nan(q_dim * AgentQuantity,1);
			for AgentNr = 1:AgentQuantity
				epsilon_k((AgentNr - 1) * q_dim + (1:q_dim)) = e_cell{AgentNr,SystemOrderNr};
			end
			epsilon((SystemOrderNr - 1) * AgentQuantity * q_dim + (1:(AgentQuantity * q_dim))) = epsilon_k;
		end
		r = nan(AgentQuantity * q_dim,1);
		for AgentNr = 1:AgentQuantity
			r_i = r_matrix(:,AgentNr);
			r((AgentNr - 1) * q_dim + (1:q_dim)) = r_i;
		end
		z = [r;epsilon];
		rho = norm((eye(AgentQuantity) - B) * ones(AgentQuantity,1) * Fl + eta_aggregated_vector);
		rho_bar = xi^(-1) * max(norm(z), chi^(-1) * vartheta_bar);
		if rho > rho_bar
			for AgentNr = 1:AgentQuantity
				if eta_aggregated_vector(AgentNr) > eta_underline_set(AgentNr)
					trigger_set(AgentNr, t_Nr) = 1;
				end
			end
		end
	end
	% Online Model Update
	for AgentNr = 1:AgentQuantity
		x_i = x_matrix(:,AgentNr);
		AgentBidirection_NeighbourSet = Bidirection_NeighbourSet{AgentNr};
		if trigger_set(AgentNr, t_Nr) == 1
			y_i = ET_MAS_GP_Leader_UnknownDynamics(x_i,q_dim,SystemOrder) + ...
				LocalGP_set{AgentNr}.SigmaN * randn(y_dim,1);
			if LocalGP_set{AgentNr}.DataQuantity >= LocalGP_set{AgentNr}.MaxDataQuantity
				LocalGP_set{AgentNr}.downdateParam(1);
			end
			LocalGP_set{AgentNr}.addPoint(x_i, y_i);
			% Update Local Prediction
			[mu_cell,var_matrix,eta_matrix] = ET_MAS_GP_Leader_LocalPrediction(AgentNr, ...
				x_matrix,LocalGP_set,beta,gamma,y_dim,mu_cell,var_matrix,eta_matrix);
			% Update Aggregation
			var_row_vector_SingleAgent = var_matrix(AgentNr,:);
			mu_row_cell_SingleAgent = mu_cell(AgentNr, :);
			[f_hat_i, eta_aggregated_i] = ET_MAS_GP_Leader_GPAggregation_SingleAgent( ...
				AgentNr, AgentBidirection_NeighbourSet, ...
				var_row_vector_SingleAgent, mu_row_cell_SingleAgent, beta, gamma);

			f_hat_matrix(:,AgentNr) = f_hat_i;
			f_true_matrix(:,AgentNr) = MAS_SingleAgent_UnknownDynamics(x_i,q_dim,SystemOrder);
			eta_aggregated_vector(AgentNr) = eta_aggregated_i;
		end

		u_cell{AgentNr} = u_cell{AgentNr} - f_hat_matrix(:,AgentNr);
	end

	u_all = ET_MAS_GP_Leader_cell2vector(u_cell,q_dim);
	% Simulation
	MultiAgent_Dynamics_FunctionHandle = @(t,x_all)ET_MAS_GP_Leader_MultiAgent_DynamicFunction(...
	t,x_all,u_all,AgentQuantity,q_dim,SystemOrder);
	[t_temp_set,x_temp_set] = ode45(MultiAgent_Dynamics_FunctionHandle, ...
		[t,t + t_step],x_all);
	%
	x_next = x_temp_set(end,:)';
	x_set(:,t_Nr + 1) = x_next;
	vartheta_all_set(:,t_Nr + 1) = x_next - s_set(:,t_Nr + 1) - kron(ones(AgentQuantity,1), x_l_set(:,t_Nr + 1));
	%
	fprintf('t = %6.4f\n',t);
end
%% State
ET_MAS_GP_Leader_Plot_SimulationState;
%% Trajectory 3D
%%
norm_vartheta_set = sqrt(sum(vartheta_all_set .^2));
figure;
semilogy(t_set,norm_vartheta_set);
%% Trigger
ET_MAS_GP_Leader_Plot_Trigger;
%% Plot Final Data Set
ET_MAS_GP_Leader_Plot_GPDataSet;


