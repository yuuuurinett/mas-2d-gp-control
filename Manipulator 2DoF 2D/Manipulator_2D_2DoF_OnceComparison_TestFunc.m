function Manipulator_2D_2DoF_OnceComparison_TestFunc( ...
	EventTriggerType,SaveFolderName,SaveFileName)
%% Set System Dimension
SystemOrder = 2;
q_dim = 2;
x_dim = q_dim * SystemOrder;
%% Set Model Paramter
m1 = 1;
m2 = 1;
L1 = 1;
L2 = 1;
g = 9.8;
%% Multi Agent System Topology
AgentQuantity = 6;
LeaderQuantity = 1;
MultiAgentSystem = Manipulator_2D_2DoF_SetMASTopology( ...
	AgentQuantity,LeaderQuantity);
%
L = MultiAgentSystem.Agent_Topology.LaplacianMatrix;
B = diag(MultiAgentSystem.Agent_Leader_Topology.ConnectionMatrix);
%%
c = 10;
lambda_set = [1;1];

%
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
% 	fprintf('The controller is stable!\n');
else
	error('The controller is not stable!');
end
% xi & chi
xi = 2 * lambda_n / min_eig_Qz * norm(Pr * (L + B));
chi = sqrt((1 + norm([t,Lambda])^2) * max_eig_Pz / min_eig_Pz ) * ...
	norm(inv(L + B));
%% Time Horizon
t_start = 0;
t_end = 50;
t_step = 0.015;
t_set = t_start:t_step:t_end;
%% Set Reference Trajectory
[xl_set,xlr_set,w] = Manipulator_2D_2DoF_LeaderDynamics(t_set,L1);
Fl = w^2 * L1 * sqrt(q_dim);
s_all_set = nan(x_dim * AgentQuantity,numel(t_set));
sr_all_set = nan(q_dim * AgentQuantity,numel(t_set));
for AgentNr = 1:AgentQuantity
	[s_all_set((AgentNr - 1) * x_dim + (1:x_dim),:), sr_all_set((AgentNr - 1) * q_dim + (1:q_dim),:)] = ...
		Manipulator_2D_2DoF_RelativePositionDynamics(t_set,AgentNr,AgentQuantity);
end
% Show Reference
% Manipulator_2D_2DoF_ShowReference(t_set,xl_set,s_all_set,AgentQuantity,L1,L2);
%% Set Gaussian Processes
% Common Setting
SigmaF = 1;
SigmaL = 0.5 * ones(x_dim,1);
GP_tau = 1e-8;
GP_delta = 0.01;
y_dim = q_dim;
LocalGP_Quantity = AgentQuantity;

DomainScale = 1.5;
X_min = DomainScale * [-1,-1,-1,-1];
X_max = DomainScale * [ 1, 1, 1, 1];
% Local_X_min = DomainScale * [ 0, 0; -1, 0; -1, -1; 0,-1]; % Local_X_min = DomainScale * [-1,-1; -1,-1; -1, -1;-1,-1];
% Local_X_max = DomainScale * [ 1, 1;  0, 1;  0,  0; 1, 0]; % Local_X_max = DomainScale * [ 1, 1;  1, 1;  1,  1; 1, 1];

% Individual Setting
MaxDataQuantity_set = 350 * ones(AgentQuantity,1);
SigmaN_set = 5 * [0.01; 0.01; 0.01; 0.01; 0.01; 0.01];

if strcmpi(EventTriggerType,'offline')
	OfflineDataQuantity_set = 1 * MaxDataQuantity_set;
else
	OfflineDataQuantity_set = 0 * MaxDataQuantity_set;
end
% if strcmpi(EventTriggerType,'distributed') || strcmpi(EventTriggerType,'centralized')
% 	OfflineDataQuantity_set = 0*[50;50;50;50];
% else
% 	OfflineDataQuantity_set = 4*[50;50;50;50];
% end

% Set Local GP
LocalGP_set = cell(LocalGP_Quantity,1);
for LocalGP_Nr = 1:LocalGP_Quantity
	MaxDataQuantity = MaxDataQuantity_set(LocalGP_Nr);
	OfflineDataQuantity = OfflineDataQuantity_set(LocalGP_Nr);
	SigmaN = SigmaN_set(LocalGP_Nr);
	LocalGP_set{LocalGP_Nr} = LocalGP_MultiOutput(x_dim,y_dim,MaxDataQuantity, ...
		SigmaN,SigmaF,SigmaL);

	X_in = 2 * (rand(x_dim,OfflineDataQuantity) - 0.5) * DomainScale;

	Y_in = Manipulator_2D_2DoF_UnknownDynamics(X_in);
	Y_in = Y_in + SigmaN * randn(size(Y_in));
	LocalGP_set{LocalGP_Nr}.add_Alldata(X_in,Y_in);

	LocalGP_set{LocalGP_Nr}.tau = GP_tau;
	LocalGP_set{LocalGP_Nr}.delta = GP_delta;
	LocalGP_set{LocalGP_Nr}.xMax = X_max;
	LocalGP_set{LocalGP_Nr}.xMin = X_min;
end

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
	warning off;
	[~,~,~,beta_i,~,~] = LocalGP_set{LocalGP_Nr}.predict(zeros(x_dim,1));
	warning on;
	beta = max(beta, beta_i);
end
beta = 1.0 * beta;
eta_underline_set = sqrt(beta) * Sigma_update_aggregation_set + gamma;
vartheta_bar = xi * chi * norm((eye(AgentQuantity) - B) * ones(AgentQuantity,1) * Fl + eta_underline_set);
%%
% rng(0);
x_all_0 = rand(x_dim * AgentQuantity,1);
x_all_set = nan(x_dim * AgentQuantity,numel(t_set));
x_all_set(:,1) = x_all_0;
vartheta_all_set = nan(AgentQuantity * x_dim,numel(t_set));
vartheta_all_set(:,1) = x_all_0 - s_all_set(:,1) - kron(ones(AgentQuantity,1), xl_set(:,1));
%
% phi_cell = cell(AgentQuantity,1);
mu_cell = cell(AgentQuantity, AgentQuantity);
var_matrix = nan(AgentQuantity, AgentQuantity);
eta_matrix = nan(AgentQuantity, AgentQuantity);
eta_aggregated_vector = nan(AgentQuantity,1);
f_hat_matrix = nan(y_dim, AgentQuantity);
f_true_matrix = nan(y_dim, AgentQuantity);
trigger_set = zeros(AgentQuantity,numel(t_set));

%%
for t_Nr = 1:numel(t_set)-1
	t = t_set(t_Nr);
% 	xl = xl_set(:,t_Nr);
	x_l_r = xlr_set(:,t_Nr);
	x_all = x_all_set(:,t_Nr);
	x_all_matrix = reshape(x_all, [], AgentQuantity);
	x_all_cell = ET_MAS_GP_Leader_vector2cell(x_all, AgentQuantity, 1);
	s_all = s_all_set(:,t_Nr);
	s_r_all = sr_all_set(:,t_Nr);
	s_r_cell = ET_MAS_GP_Leader_vector2cell(s_r_all, AgentQuantity, 1);
	x_tilde_all = x_all - s_all;
	x_tilde_cell = ET_MAS_GP_Leader_vector2cell(x_tilde_all, AgentQuantity, SystemOrder);
	vartheta_all = vartheta_all_set(:,t_Nr);
	vartheta_cell = ET_MAS_GP_Leader_vector2cell(vartheta_all, AgentQuantity, SystemOrder);
	% Consensus Control
	[phi_cell,r_matrix,e_cell] = Manipulator_2D_2DoF_ConsensusLaw(vartheta_cell,x_tilde_cell,x_l_r, ...
		MultiAgentSystem,c,lambda_set,s_r_cell);
	% Learning Part
	if strcmpi(EventTriggerType,'exact')
	for AgentNr = 1:AgentQuantity
		x_i = x_all_cell{AgentNr};
		[~,~,f_i] = Manipulator_2D_2DoF_get_Dynamics_h_g_f(x_i,L1,L2,m1,m2);
		f_hat_matrix(:,AgentNr) = f_i;
	end
	else
	for AgentNr = 1:AgentQuantity
		x_i = x_all_matrix(:,AgentNr);
		% Local prediction
		[mu_cell{AgentNr, AgentNr},var_matrix(AgentNr, AgentNr),eta_matrix(AgentNr, AgentNr)] = ...
			Manipulator_2D_2DoF_LocalPrediction(x_i,AgentNr,LocalGP_set,beta,gamma,y_dim);
		% Cooperative Prediction
		AgentBidirection_NeighbourSet = Bidirection_NeighbourSet{AgentNr};
		AgentBidirection_NeighbourQuantity = numel(AgentBidirection_NeighbourSet);
		for Bidirection_NeighborNr = 1:AgentBidirection_NeighbourQuantity
			Bidirection_NeighborAgentNr = AgentBidirection_NeighbourSet(Bidirection_NeighborNr);
			[mu_cell{AgentNr, Bidirection_NeighborAgentNr}, ...
				var_matrix(AgentNr, Bidirection_NeighborAgentNr), ...
				eta_matrix(AgentNr, Bidirection_NeighborAgentNr)] = ...
				Manipulator_2D_2DoF_LocalPrediction( ...
				x_i,Bidirection_NeighborAgentNr,LocalGP_set,beta,gamma,y_dim);
		end
		% GP aggregation - Posterior variance
		var_row_vector_SingleAgent = var_matrix(AgentNr,:);
		mu_row_cell_SingleAgent = mu_cell(AgentNr, :);
		[f_hat_i, eta_aggregated_i] = ET_MAS_GP_Leader_GPAggregation_SingleAgent( ...
			AgentNr, AgentBidirection_NeighbourSet, ...
			var_row_vector_SingleAgent, mu_row_cell_SingleAgent, beta, gamma);
		% Save Distributed GP result
		f_hat_matrix(:,AgentNr) = f_hat_i;
		f_true_matrix(:,AgentNr) = Manipulator_2D_2DoF_UnknownDynamics(x_i);
		eta_aggregated_vector(AgentNr) = eta_aggregated_i;
	end
	%% Event-trigger
	for AgentNr = 1:AgentQuantity
		switch EventTriggerType
			case 'distributed'
				trigger_set(AgentNr, t_Nr) = ...
					Manipulator_2D_2DoF_DistributedET(AgentNr, ...
					r_matrix,e_cell,eta_underline_set,eta_aggregated_vector, ...
					MultiAgentSystem,Fl,xi,chi,vartheta_bar);
			case 'centralized'
				trigger_flag = Manipulator_2D_2DoF_CentralizedET( ...
					r_matrix,e_cell,eta_aggregated_vector, ...
					MultiAgentSystem,Fl,xi,chi,vartheta_bar);
				if trigger_flag == 1 && eta_aggregated_vector(AgentNr) > eta_underline_set(AgentNr)
					trigger_set(AgentNr, t_Nr) = 1;
				end
			case 'time'
				trigger_set(AgentNr, t_Nr) = 1;
			case 'offline'
				trigger_set(AgentNr, t_Nr) = 0;
		end
	end
	%% Online Learning
	for AgentNr = 1:AgentQuantity
		x_i = x_all_matrix(:,AgentNr);
		AgentBidirection_NeighbourSet = Bidirection_NeighbourSet{AgentNr};
		if trigger_set(AgentNr, t_Nr) == 1
			y_i = Manipulator_2D_2DoF_UnknownDynamics(x_i) + ...
				LocalGP_set{AgentNr}.SigmaN * randn(y_dim,1);
			if LocalGP_set{AgentNr}.DataQuantity >= LocalGP_set{AgentNr}.MaxDataQuantity
				LocalGP_set{AgentNr}.downdateParam(1);
			end
			LocalGP_set{AgentNr}.addPoint(x_i, y_i);
			% Update Local Prediction
			[mu_cell{AgentNr, AgentNr},var_matrix(AgentNr, AgentNr),eta_matrix(AgentNr, AgentNr)] = ...
				Manipulator_2D_2DoF_LocalPrediction(x_i,AgentNr,LocalGP_set,beta,gamma,y_dim);
			% Update Aggregation
			var_row_vector_SingleAgent = var_matrix(AgentNr,:);
			mu_row_cell_SingleAgent = mu_cell(AgentNr, :);
			[f_hat_i, eta_aggregated_i] = ET_MAS_GP_Leader_GPAggregation_SingleAgent( ...
				AgentNr, AgentBidirection_NeighbourSet, ...
				var_row_vector_SingleAgent, mu_row_cell_SingleAgent, beta, gamma);

			f_hat_matrix(:,AgentNr) = f_hat_i;
			eta_aggregated_vector(AgentNr) = eta_aggregated_i;
		end
	end
	end
	%% Calculate Control Input
	u_cell = Manipulator_2D_2DoF_get_u_cell(x_all_cell,phi_cell,f_hat_matrix,L1,L2,m1,m2);
% 	u_all = ET_MAS_GP_Leader_cell2vector(u_cell,q_dim);
	%% Simulation
	[~,x_all_temp] = ode45( ...
		@(t,x_all)Manipulator_2D_2DoF_MultiAgent_DynamicFunction(t,x_all,u_cell,L1,L2,m1,m2), ...
		[t, t + t_step],x_all);
	x_all_next = x_all_temp(end,:)';
	x_all_set(:,t_Nr + 1) = x_all_next;
	vartheta_all_set(:,t_Nr + 1) = x_all_next - s_all_set(:,t_Nr + 1) - ...
		kron(ones(AgentQuantity,1), xl_set(:,t_Nr + 1));
	%% Print Time
% 	fprintf('t = %6.4f\n',t);
end
save([SaveFolderName,'\',SaveFileName,'.mat']);

% save(['Result\Event Trigger for MAS using GP with Leader\Manipulator 2DoF 2D\Once Comparison\', ...
% 	EventTriggerType,'.mat']);
end